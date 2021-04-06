//
// Created by per on 03.02.21.
//

#include "sortledton_driver.hpp"

#include <chrono>
#include <assert.h>
#include <fstream>
#include <omp.h>

#include "common/timer.hpp"
#include "third-party/gapbs/gapbs.hpp"
#include "third-party/libcuckoo/cuckoohash_map.hh"
#include "utility/timeout_service.hpp"

#include "not_implemented.hpp"

#include "data_types.h"
#include "versioning/EdgeDoesNotExistsPrecondition.h"
#include "versioning/VersionedEdgeIterator.h"

using namespace gapbs;

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG
namespace gfe { extern mutex _log_mutex [[maybe_unused]]; }
#define COUT_DEBUG_FORCE(msg) { std::scoped_lock<std::mutex> lock{::gfe::_log_mutex}; std::cout << "[TeseoDriver::" << __FUNCTION__ << "] " << msg << std::endl; }
#if defined(DEBUG)
#define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
#define COUT_DEBUG(msg)
#endif


namespace gfe::library {

    SortledtonDriver::SortledtonDriver(bool is_graph_directed, bool sparse_graph, uint64_t max_num_vertices,
                                       int block_size) : tm(1), m_is_directed(is_graph_directed) {
      if (is_graph_directed == true) {
        throw std::invalid_argument("Only undirected graphs are currently supported by the front-end");
      }
      if (sparse_graph == true) { throw std::invalid_argument("Only dense graphs are supported by the front-end."); }
      ds = new VersioningBlockedSkipListAdjacencyList(block_size, tm);
      ds->reserve_vertices(max_num_vertices);
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
      SnapshotTransaction tx = non_const_this->tm.getSnapshotTransaction(ds);
      auto num_edges = tx.edge_count() / 2;
      non_const_this->tm.transactionCompleted(tx);
      return num_edges;
    }

    uint64_t SortledtonDriver::num_vertices() const {
      SortledtonDriver *non_const_this = const_cast<SortledtonDriver *>(this);
      SnapshotTransaction tx = non_const_this->tm.getSnapshotTransaction(ds);
      auto num_vertices = tx.vertex_count();
      non_const_this->tm.transactionCompleted(tx);
      return num_vertices;
    }

/**
 * Returns true if the given vertex is present, false otherwise
 */
    bool SortledtonDriver::has_vertex(uint64_t vertex_id) const {
      SortledtonDriver *non_const_this = const_cast<SortledtonDriver *>(this);
      SnapshotTransaction tx = non_const_this->tm.getSnapshotTransaction(ds);  // TODO weights currently not supported
      auto has_vertex = tx.has_vertex(vertex_id);
      non_const_this->tm.transactionCompleted(tx);
      return has_vertex;
    }

/**
 * Returns the weight of the given edge is the edge is present, or NaN otherwise
 */
    double SortledtonDriver::get_weight(uint64_t source, uint64_t destination) const {
      SortledtonDriver *non_const_this = const_cast<SortledtonDriver *>(this);
      SnapshotTransaction tx = non_const_this->tm.getSnapshotTransaction(ds);  // TODO weights currently not supported
      auto has_edge = tx.has_edge({static_cast<dst_t>(source), static_cast<dst_t>(destination)});
      non_const_this->tm.transactionCompleted(tx);
      return has_edge ? 0.0 : nan("");
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
      SnapshotTransaction tx = tm.getSnapshotTransaction(ds);
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
      SnapshotTransaction tx = tm.getSnapshotTransaction(ds);

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
        tx.insert_edge(internal_edge);
        tx.insert_edge({internal_edge.dst, internal_edge.src});
        inserted &= tx.execute();
      } catch (exception &e) {
        inserted = false;
      }
      tm.transactionCompleted(tx);
      return inserted;
    }

    bool SortledtonDriver::add_edge_v2(gfe::graph::WeightedEdge e) {
      assert(!m_is_directed);
      edge_t internal_edge{static_cast<dst_t>(e.source()), static_cast<dst_t>(e.destination())};

      SnapshotTransaction tx = tm.getSnapshotTransaction(ds);
      tx.use_vertex_does_not_exists_semantics();

      tx.insert_vertex(internal_edge.src);
      tx.insert_vertex(internal_edge.dst);

      tx.insert_edge({internal_edge.dst, internal_edge.src});
      tx.insert_edge(internal_edge);

      bool inserted = true;
      inserted &= tx.execute();

      tm.transactionCompleted(tx);
      return inserted;
    }

    bool SortledtonDriver::remove_edge(gfe::graph::Edge e) {
      throw NotImplemented();
    }


    /*****************************************************************************
    *                                                                           *
    *  BFS                                                                      *
    *                                                                           *
    ****************************************************************************/
    namespace { // anonymous

        /*
        GAP Benchmark Suite
        Kernel: Breadth-First Search (BFS)
        Author: Scott Beamer
        Will return parent array for a BFS traversal from a source vertex
        This BFS implementation makes use of the Direction-Optimizing approach [1].
        It uses the alpha and beta parameters to determine whether to switch search
        directions. For representing the frontier, it uses a SlidingQueue for the
        top-down approach and a Bitmap for the bottom-up approach. To reduce
        false-sharing for the top-down approach, thread-local QueueBuffer's are used.
        To save time computing the number of edges exiting the frontier, this
        implementation precomputes the degrees in bulk at the beginning by storing
        them in parent array as negative numbers. Thus the encoding of parent is:
          parent[x] < 0 implies x is unvisited and parent[x] = -out_degree(x)
          parent[x] >= 0 implies x been visited
        [1] Scott Beamer, Krste Asanović, and David Patterson. "Direction-Optimizing
            Breadth-First Search." International Conference on High Performance
            Computing, Networking, Storage and Analysis (SC), Salt Lake City, Utah,
            November 2012.
        */

        static int64_t BUStep(VersioningBlockedSkipListAdjacencyList *ds, SnapshotTransaction &tx, pvector <int64_t> &distances, int64_t distance, Bitmap &front, Bitmap &next) {
          const int64_t N = tx.vertex_count();
          int64_t awake_count = 0;
          next.reset();
#pragma omp parallel for reduction(+ : awake_count) schedule(dynamic, 1024)
          for (int64_t u = 0; u < N; u++) {
            if (distances[u] < 0) { // the node has not been visited yet
              bool done = false;
              sortledton_iterator iter(*ds);
              tx.neighbourhood(u, iter);
              while (!done && iter.has_next()) {
                dst_t n = iter.next();
                if (front.get_bit(n)) {
                  distances[u] = distance; // on each BUStep, all nodes will have the same distance
                  awake_count++;
                  next.set_bit(u);
                  done = true;
                  iter.close();
                }
              }
            }
          }
          return awake_count;
        }

        static int64_t TDStep(VersioningBlockedSkipListAdjacencyList *ds, SnapshotTransaction &tx, pvector <int64_t> &distances, int64_t distance, SlidingQueue <int64_t> &queue) {
          int64_t scout_count = 0;

          #pragma omp parallel
          {
            QueueBuffer <int64_t> lqueue(queue);
            #pragma omp for reduction(+ : scout_count)
            for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
              int64_t u = *q_iter;

              sortledton_iterator iter(*ds);
              tx.neighbourhood(u, iter);
              while (iter.has_next()) {
                dst_t destination = iter.next();
                int64_t curr_val = distances[destination];

                if (curr_val < 0 && compare_and_swap(distances[destination], curr_val, distance)) {
                  lqueue.push_back(destination);
                  scout_count += -curr_val;
                }
              }
            }
            lqueue.flush();
          }
          return scout_count;
        }

        static void QueueToBitmap(const SlidingQueue <int64_t> &queue, Bitmap &bm) {
          #pragma omp parallel for
          for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
            int64_t u = *q_iter;
            bm.set_bit_atomic(u);
          }
        }

        static void BitmapToQueue(SnapshotTransaction &tx, const Bitmap &bm, SlidingQueue <int64_t> &queue) {
          const int64_t N = tx.vertex_count();

          #pragma omp parallel
          {
            QueueBuffer <int64_t> lqueue(queue);
            #pragma omp for
            for (int64_t n = 0; n < N; n++) {
              if (bm.get_bit(n)) {
                lqueue.push_back(n);
              }
            }
            lqueue.flush();
          }
          queue.slide_window();
        }

        static pvector <int64_t> InitDistances(SnapshotTransaction &tx) {
          const int64_t N = tx.vertex_count();
          pvector <int64_t> distances(N);

          #pragma omp parallel for
          for (int64_t n = 0; n < N; n++) {
            int64_t out_degree = tx.neighbourhood_size(n);
            distances[n] = out_degree != 0 ? -out_degree : -1;
          }
          return distances;
        }

    } // anon namespace

    static pvector <int64_t> sortledton_bfs(VersioningBlockedSkipListAdjacencyList* ds, SnapshotTransaction& tx, int64_t source, utility::TimeoutService &timer, int alpha = 15, int beta = 18) {
      // The implementation from GAP BS reports the parent (which indeed it should make more sense), while the one required by
      // Graphalytics only returns the distance

      pvector <int64_t> distances = InitDistances(tx);
      distances[source] = 0;

      uint64_t vertex_count = tx.vertex_count();
      SlidingQueue <int64_t> queue(vertex_count);
      queue.push_back(source);
      queue.slide_window();
      Bitmap curr(vertex_count);
      curr.reset();
      Bitmap front(vertex_count);
      front.reset();
      int64_t edges_to_check = tx.edge_count();  // TODO this could be a slow down given my implementation, we could sum up the adjacency set sizes in a parallel for loop
      int64_t scout_count = tx.neighbourhood_size(source);
      int64_t distance = 1; // current distance
      while (!timer.is_timeout() && !queue.empty()) {
        if (scout_count > edges_to_check / alpha) {
          int64_t awake_count, old_awake_count;
          QueueToBitmap(queue, front);
          awake_count = queue.size();
          queue.slide_window();
          do {
            old_awake_count = awake_count;
            awake_count = BUStep(ds, tx, distances, distance, front, curr);
            front.swap(curr);
            distance++;
          } while ((awake_count >= old_awake_count) ||
                   (awake_count > static_cast<int64_t>(vertex_count) / beta));
          BitmapToQueue(tx, front, queue);
          scout_count = 1;
        } else {
          edges_to_check -= scout_count;
          scout_count = TDStep(ds, tx, distances, distance, queue);
          queue.slide_window();
          distance++;
        }
      }

      return distances;
    }

    static void save_bfs(libcuckoo::cuckoohash_map <uint64_t, int64_t> &result, const char *dump2file) {
      assert(dump2file != nullptr);
      COUT_DEBUG("save the results to: " << dump2file)

      fstream handle(dump2file, ios_base::out);
      if (!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

      auto list_entries = result.lock_table();

      for (const auto &p : list_entries) {
        handle << p.first << " ";

        // if  the vertex was not reached, the algorithm sets its distance to < 0
        if (p.second < 0) {
          handle << numeric_limits<int64_t>::max();
        } else {
          handle << p.second;
        }
        handle << "\n";
      }

      list_entries.unlock();
      handle.close();
    }

    static void check_bfs(vertex_id_t start_vertex, vector <uint> &distances) {
      const string gold_standard_directory = "/space/fuchs/shared/graph_two_gold_standards";

      const string gold_standard_file =
              gold_standard_directory + "/bfs_" + "live-journal-full-insert" + "_" + "0" + "_" +
              "inserts" + ".goldStandard";

      cout << "Validating bfs experiment against " << gold_standard_file << endl;
      ifstream f(gold_standard_file, ifstream::in | ifstream::binary);

      size_t size;
      f.read((char *) &size, sizeof(size));
      cout << "Expected size " << size << endl;
      cout << "Actual size " << distances.size();
//      assert(size == distances.size());

      uint e;
      int i = 0;
      int errors = 0;
      for (auto d : distances) {
        f.read((char *) &e, sizeof(e));
        if (d != e) {
          errors++;
//          if (errors < 100) {
            cout << "i " << i << " d " << d << " e " << e << endl;
//          }
        }
//        assert(d == e);
        i += 1;
      }
      cout << "Total number of errors " << errors << endl;

      f.close();
    }


    bool SortledtonDriver::gced = false;

    void SortledtonDriver::bfs(uint64_t source_vertex_id, const char *dump2file) {
      utility::TimeoutService tcheck{m_timeout};
      common::Timer timer;
      timer.start();

      // TODO GC could be parallelized
      if (!gced) {
        ds->gc_all();
        cout << "gc done after " << timer << endl;
        gced = true;
      }
      SnapshotTransaction tx = tm.getSnapshotTransaction(ds);

      // execute the BFS algorithm
      auto distances = sortledton_bfs(ds, tx, source_vertex_id, tcheck);
      if (tcheck.is_timeout()) { RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }
      cout << "BFS took " << timer << endl;

      int N = distances.size();

      // Check bfs algorithm
//      vector <uint> ret(N);
//#pragma omp parallel for
//      for (int i = 0; i < N; i++) {
//        if (distances[i] < 0) {
//          ret[i] = numeric_limits<uint>::max();
//        } else {
//          ret[i] = distances[i];
//        }
//      }
//      check_bfs(0, ret);

      libcuckoo::cuckoohash_map </* external id */ uint64_t, /* distance */ int64_t> external_ids;
      #pragma omp parallel for
      for (uint64_t i = 0; i < N; i++) {
        uint64_t external_node_id = i;

        auto distance = distances[i];

//        cout << external_node_id << " " << distance << endl;

        external_ids.insert(external_node_id, distance);
      }
      tm.transactionCompleted(tx); // TODO do this in destructor finally.

      if (tcheck.is_timeout()) RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);

      // store the results in the given file
      if (dump2file != nullptr) {
        save_bfs(external_ids, dump2file); // TODO what needs to be saved
      }
    }



    void SortledtonDriver::pagerank(uint64_t num_iterations, double damping_factor, const char *dump2file) {
      throw NotImplemented();
    }

    void SortledtonDriver::wcc(const char *dump2file) {
      throw NotImplemented();
    }

    void SortledtonDriver::cdlp(uint64_t max_iterations, const char *dump2file) {
      throw NotImplemented();
    }

    void SortledtonDriver::lcc(const char *dump2file) {
      throw NotImplemented();
    }

    void SortledtonDriver::sssp(uint64_t source_vertex_id, const char *dump2file) {
      throw NotImplemented();
    }

    bool SortledtonDriver::can_be_validated() const {
     return true;
    }
}