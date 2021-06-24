//
// Created by per on 03.02.21.
//

#include "sortledton_driver.hpp"

#include <chrono>
#include <optional>
#include <omp.h>

#include "third-party/gapbs/gapbs.hpp"

#include "common/timer.hpp"
#include "third-party/gapbs/gapbs.hpp"
#include "utility/timeout_service.hpp"

#include "not_implemented.hpp"

#include "data_types.h"
#include "versioning/EdgeDoesNotExistsPrecondition.h"
#include "versioning/VersionedEdgeIterator.h"
#include "experiments/SSSP.h"
#include "experiments/PageRank.h"
#include "experiments/WCC.h"
#include "experiments/CDLP.h"
#include "experiments/LCC.h"
#include "experiments/GAPBSAlgorithms.h"

using namespace gapbs;
using namespace common;

namespace gfe { extern mutex _log_mutex [[maybe_unused]]; }


namespace gfe::library {

    SortledtonDriver::SortledtonDriver(bool is_graph_directed, size_t properties_size, int block_size) : tm(1), m_is_directed(is_graph_directed) {
      if (is_graph_directed == true) {
        throw std::invalid_argument("Only undirected graphs are currently supported by the front-end");
      }
      ds = new VersioningBlockedSkipListAdjacencyList(block_size, properties_size, tm);
    }

    SortledtonDriver::~SortledtonDriver() {
      delete ds;
      ds = nullptr;
    }

    void SortledtonDriver::on_main_init(int num_threads) {
      tm.reset_max_threads(num_threads);
    }

    void SortledtonDriver::on_thread_init(int thread_id) {
      tm.register_thread(thread_id);
    }

    void SortledtonDriver::on_thread_destroy(int thread_id) {
      tm.deregister_thread(thread_id);
    }

    void SortledtonDriver::dump_ostream(std::ostream &out) const {
      throw NotImplemented();
    }

    uint64_t SortledtonDriver::num_edges() const {
      SortledtonDriver *non_const_this = const_cast<SortledtonDriver *>(this);
      SnapshotTransaction tx = non_const_this->tm.getSnapshotTransaction(ds, false);
      auto num_edges = tx.edge_count() / 2;
      non_const_this->tm.transactionCompleted(tx);
      return num_edges;
    }

    uint64_t SortledtonDriver::num_vertices() const {
      SortledtonDriver *non_const_this = const_cast<SortledtonDriver *>(this);
      SnapshotTransaction tx = non_const_this->tm.getSnapshotTransaction(ds, false);
      auto num_vertices = tx.vertex_count();
      non_const_this->tm.transactionCompleted(tx);
      return num_vertices;
    }

/**
 * Returns true if the given vertex is present, false otherwise
 */
    bool SortledtonDriver::has_vertex(uint64_t vertex_id) const {
      SortledtonDriver *non_const_this = const_cast<SortledtonDriver *>(this);
      SnapshotTransaction tx = non_const_this->tm.getSnapshotTransaction(ds, false);  // TODO weights currently not supported
      auto has_vertex = tx.has_vertex(vertex_id);
      non_const_this->tm.transactionCompleted(tx);
      return has_vertex;
    }

/**
 * Returns the weight of the given edge is the edge is present, or NaN otherwise
 */
    double SortledtonDriver::get_weight(uint64_t source, uint64_t destination) const {
      SortledtonDriver *non_const_this = const_cast<SortledtonDriver *>(this);
      SnapshotTransaction tx = non_const_this->tm.getSnapshotTransaction(ds, false);  // TODO weights currently not supported
      weight_t w;
      auto has_edge = tx.get_weight({static_cast<dst_t>(source), static_cast<dst_t>(destination)}, (char*) &w);
      non_const_this->tm.transactionCompleted(tx);
      return has_edge ? w : nan("");
    }

/**
 * Check whether the graph is directed
 */
    bool SortledtonDriver::is_directed() const {
      return m_is_directed;
    }

/**
 * Impose a timeout on each graph computation. A computation that does not terminate by the given seconds will raise a TimeoutError.
 */
    void SortledtonDriver::set_timeout(uint64_t seconds) {
      m_timeout = chrono::seconds{seconds};
    }

/**
 * Add the given vertex to the graph
 * @return true if the vertex has been inserted, false otherwise (that is, the vertex already exists)
 */
    bool SortledtonDriver::add_vertex(uint64_t vertex_id) {
      SnapshotTransaction tx = tm.getSnapshotTransaction(ds, true);
      bool inserted = true;
      try {
        tx.insert_vertex(vertex_id);
        tx.execute();
      } catch (exception &e) {
        inserted = false;
      }
      tm.transactionCompleted(tx);
      return inserted;
    }

/**
 * Remove the mapping for a given vertex. The actual internal vertex is not removed from the adjacency list.
 * @param vertex_id the vertex to remove
 * @return true if a mapping for that vertex existed, false otherwise
 */
    bool SortledtonDriver::remove_vertex(uint64_t vertex_id) {
      throw NotImplemented();
    }


/**
 * Adds a given edge to the graph if both vertices exists already
 */
    bool SortledtonDriver::add_edge(gfe::graph::WeightedEdge e) {
      assert(!m_is_directed);
      edge_t internal_edge{static_cast<dst_t>(e.source()), static_cast<dst_t>(e.destination())};

      thread_local optional<SnapshotTransaction> tx_o = nullopt;
      if (tx_o.has_value()) {
        tm.getSnapshotTransaction(ds, true, *tx_o);
      } else {
        tx_o = tm.getSnapshotTransaction(ds, true);
      }
      auto tx = *tx_o;

      VertexExistsPrecondition pre_v1(internal_edge.src);
      tx.register_precondition(&pre_v1);
      VertexExistsPrecondition pre_v2(internal_edge.dst);
      tx.register_precondition(&pre_v2);
      // Even in the undirected case, we need to check only for the existence of one edge direction to ensure consistency.
      EdgeDoesNotExistsPrecondition pre_e(internal_edge);
      tx.register_precondition(&pre_e);

      // test
      bool inserted = true;
      try {
        tx.insert_edge(internal_edge, (char*) &e.m_weight, sizeof(e.m_weight));
        tx.insert_edge({internal_edge.dst, internal_edge.src}, (char*) &e.m_weight, sizeof(e.m_weight));
        inserted &= tx.execute();
      } catch (exception &e) {
        inserted = false;
      }
      tm.transactionCompleted(tx);
      return inserted;
    }

    bool SortledtonDriver::add_edge_v2(gfe::graph::WeightedEdge e) {
      assert(!m_is_directed);

      thread_local optional<SnapshotTransaction> tx_o = nullopt;
      if (tx_o.has_value()) {
        tm.getSnapshotTransaction(ds, true, *tx_o);
      } else {
        tx_o = tm.getSnapshotTransaction(ds, true);
      }
      auto tx = *tx_o;

      edge_t internal_edge{static_cast<dst_t>(e.source()), static_cast<dst_t>(e.destination())};

      tx.use_vertex_does_not_exists_semantics();

      tx.insert_vertex(internal_edge.src);
      tx.insert_vertex(internal_edge.dst);

      tx.insert_edge({internal_edge.dst, internal_edge.src}, (char*) &e.m_weight, sizeof(e.m_weight));
      tx.insert_edge(internal_edge, (char*) &e.m_weight, sizeof(e.m_weight));

      bool inserted = true;
      inserted &= tx.execute();

      tm.transactionCompleted(tx);
      return inserted;
    }

    bool SortledtonDriver::remove_edge(gfe::graph::Edge e) {
      assert(!m_is_directed);

      thread_local optional<SnapshotTransaction> tx_o = nullopt;
      if (tx_o.has_value()) {
        tm.getSnapshotTransaction(ds, true, *tx_o);
      } else {
        tx_o = tm.getSnapshotTransaction(ds, true);
      }
      auto tx = *tx_o;

      edge_t internal_edge{static_cast<dst_t>(e.source()), static_cast<dst_t>(e.destination())};

      tx.delete_edge({internal_edge.dst, internal_edge.src});
      tx.delete_edge(internal_edge);

      bool removed = true;
      removed &= tx.execute();

      tm.transactionCompleted(tx);
      return removed;
      ;
    }

    void SortledtonDriver::run_gc() {
      if (!gced) {
        Timer t; t.start();
        ds->gc_all();
        gced = true;
        cout << "Running GC took: " << t;
      }
    }

    static void save_bfs(vector<pair<uint64_t, uint>>& result, const char *dump2file) {
      assert(dump2file != nullptr);
      COUT_DEBUG("save the results to: " << dump2file)

      fstream handle(dump2file, ios_base::out);
      if (!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

      for (const auto &p : result) {
        handle << p.first << " ";

        // if  the vertex was not reached, the algorithm sets its distance to < 0
        if (p.second == numeric_limits<uint>::max()) {
          handle << numeric_limits<int64_t>::max();
        } else {
          handle << (int64_t) p.second;
        }
        handle << "\n";
      }
      handle.close();
    }

    static vector<pair<uint64_t, uint>> translate_bfs(SnapshotTransaction& tx, pvector<int64_t>& values) {
      auto N = values.size();

      vector<pair<vertex_id_t , uint>> logical_result(N);

#pragma omp parallel for
      for (uint v = 0; v <  N; v++) {
        if (tx.has_vertex_p(v)) {
          if (values[v] >= 0) {
            logical_result[v] = make_pair(tx.logical_id(v), values[v]);
          } else {
            logical_result[v] = make_pair(tx.logical_id(v), numeric_limits<uint>::max());
          }
        } else {
          logical_result[v] = make_pair(v, numeric_limits<uint>::max());
        }
      }
      return logical_result;
    }


    void SortledtonDriver::bfs(uint64_t source_vertex_id, const char *dump2file) {
      tm.register_thread(0);
      run_gc();

      SnapshotTransaction tx = tm.getSnapshotTransaction(ds, false);
      auto physical_src = tx.physical_id(source_vertex_id);

      Timer t; t.start();
      auto distances = GAPBSAlgorithms::bfs(tx, physical_src, false);

      cout << "BFS took " << t << endl;
      auto external_ids = translate_bfs(tx, distances);
      cout << "Translation took " << t << endl;
      tm.transactionCompleted(tx);

      if (dump2file != nullptr) {
        save_bfs(external_ids, dump2file);
      }
      tm.deregister_thread(0);
    }



    void SortledtonDriver::pagerank(uint64_t num_iterations, double damping_factor, const char *dump2file) {
      tm.register_thread(0);
      run_gc();

      SnapshotTransaction tx = tm.getSnapshotTransaction(ds, false);

      auto pr = PageRank::page_rank_bs(tx, num_iterations, damping_factor);;
      auto external_ids = translate<double>(tx, pr);

      tm.transactionCompleted(tx);

      if (dump2file != nullptr) {
        save_result<double>(external_ids, dump2file);
      }
      tm.deregister_thread(0);
    }

    void SortledtonDriver::wcc(const char *dump2file) {
      tm.register_thread(0);
      run_gc();

      SnapshotTransaction tx = tm.getSnapshotTransaction(ds, false);

      auto clusters = WCC::gapbs_wcc(tx);
      auto external_ids = translate<uint64_t>(tx, clusters);

      tm.transactionCompleted(tx);

      if (dump2file != nullptr) {
        save_result<uint64_t>(external_ids, dump2file);
      }
      tm.deregister_thread(0);
    }

    void SortledtonDriver::cdlp(uint64_t max_iterations, const char *dump2file) {
      tm.register_thread(0);
      run_gc();

      SnapshotTransaction tx = tm.getSnapshotTransaction(ds, false);

      auto clusters = CDLP::teseo_cdlp(tx, max_iterations);
      auto external_ids = translate<uint64_t>(tx, clusters);

      tm.transactionCompleted(tx);

      if (dump2file != nullptr) {
        save_result<uint64_t>(external_ids, dump2file);
      }
      tm.deregister_thread(0);
    }

    void SortledtonDriver::lcc(const char *dump2file) {
      tm.register_thread(0);
      run_gc();

      SnapshotTransaction tx = tm.getSnapshotTransaction(ds, false);

      auto lcc_values = LCC::lcc_merge_sort(tx);
      auto external_ids = translate<double>(tx, lcc_values);

      tm.transactionCompleted(tx);

      if (dump2file != nullptr) {
        save_result<double>(external_ids, dump2file);
      }
      tm.deregister_thread(0);
    }

    void SortledtonDriver::sssp(uint64_t source_vertex_id, const char *dump2file) {
      tm.register_thread(0);
      run_gc();

      SnapshotTransaction tx = tm.getSnapshotTransaction(ds, false);
      auto physical_src = tx.physical_id(source_vertex_id);

      auto distances = SSSP::gabbs_sssp(tx, physical_src, 2.0);

      auto external_ids = translate<double>(tx, distances);

      tm.transactionCompleted(tx);

      if (dump2file != nullptr) {
        save_result<double>(external_ids, dump2file);
      }
      tm.deregister_thread(0);
    }

    bool SortledtonDriver::can_be_validated() const {
     return true;
    }
}