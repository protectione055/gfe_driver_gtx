//
// Created by Ziming Zhang on 8/13/24.
//

#include "age_driver.hpp"

#include <atomic>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <mutex>
#include <limits>
#include <sstream>
#include <thread>
#include <unordered_set>
#include <stdexcept>
#include <algorithm>

#include "../../third-party/libcommon/include/lib/common/system.hpp"
#include "../../third-party/libcommon/include/lib/common/timer.hpp"
#include "../../third-party/gapbs/gapbs.hpp"
#include "../../utility/timeout_service.hpp"

using namespace common;
using namespace std;

using vertex_t = uint64_t;

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG
namespace gfe { extern mutex _log_mutex [[maybe_unused]]; }

#define LOG(level, msg) { std::scoped_lock<std::mutex> lock{::gfe::_log_mutex}; std::cerr << "["#level"]" << "[AGEDriver::" << __FUNCTION__ << "] [Thread #" << common::concurrency::get_thread_id() << "] " << msg << std::endl; }

#define GRAPH_NAME this->_graph_name

/*****************************************************************************
 *                                                                           *
 *  Common Utilities                                                         *
 *                                                                           *
 *****************************************************************************/
#define DBL_EQ(v1, v2) (fabs(v1 - v2) < DBL_EPSILON)

namespace gfe::library {
    // ConnectionPool
    ConnectionPool::ConnectionPool(const std::string& conninfo, std::size_t pool_size)
     : conninfo_(conninfo), pool_size_(pool_size) {
        for (std::size_t i = 0; i < pool_size_; ++i) {
            try {
                    connections_.push(create_connection(conninfo));
            } catch (const std::exception& e) {
                std::cerr << e.what() << std::endl;
            }
        }
    }

    ConnectionPool::~ConnectionPool() {
        while (!connections_.empty()) {
            auto conn = connections_.front();
            connections_.pop();
            conn->disconnect();
        }
    }

    std::unique_ptr<pqxx::connection> ConnectionPool::create_connection(const std::string& conninfo) {
        std::unique_ptr<pqxx::connection> conn = std::make_unique<pqxx::connection>(conninfo);
        if (conn->is_open()) {
            pqxx::work txn(*conn);

            // LOAD plugin
            pqxx::result res = txn.exec("LOAD 'age'");

            // SET search_path
            res = txn.exec( "SELECT pg_catalog.set_config('search_path', 'ag_catalog, \"$user\", public', false)");

            txn.commit();
        }
        return conn;
    }

    std::shared_ptr<pqxx::connection> ConnectionPool::acquire() {
        std::unique_lock<std::mutex> lock(mutex_);
        condition_.wait(lock, [this]() { return !connections_.empty(); });
        
        auto conn = connections_.front();
        connections_.pop();
        return conn;
    }

    void ConnectionPool::release(std::shared_ptr<pqxx::connection> conn) {
        std::unique_lock<std::mutex> lock(mutex_);
        connections_.push(conn);
        condition_.notify_one();
    }

    void ConnectionPool::resize(size_t size) {
        std::unique_lock<std::mutex> lock(mutex_);
        pool_size_ = size;
        while (connections_.size() > pool_size_) {
            auto conn = connections_.front();
            connections_.pop();
            release(conn);
        }
        while (connections_.size() < pool_size_) {
            try {
                    connections_.push(create_connection(conninfo_));
            } catch (const std::exception& e) {
                LOG(ERROR, e.what());
                exit(1);
            }
        }
    }

    // AGEDriver
    AGEDriver::AGEDriver(bool is_directed, bool read_only): m_is_directed(is_directed),m_read_only(read_only) { 
        connect_to_database(8, "localhost", "5432", "postgres", "zzm", "66668888");
        open_graph("gfe_driver");
    }

    AGEDriver::~AGEDriver() noexcept {
        // Do nothing
    }

    void AGEDriver::connect_to_database(size_t pool_size, const std::string& host, const std::string& port, const std::string& dbname, const std::string& user, const std::string& password) {
        std::ostringstream conninfo;
        conninfo << "host" << "=" << host << " ";
        conninfo << "port" << "=" << port << " ";
        conninfo << "dbname" << "=" << dbname << " ";
        conninfo << "user" << "=" << user << " ";
        conninfo << "password" << "=" << password << " ";

        _conn_pool = std::make_unique<ConnectionPool>(conninfo.str(), pool_size);
    }

    void AGEDriver::open_graph(const std::string& graph_name) {
        /*
         * Check if graph_name already exists
         * If not exists, 
         */
        auto conn = _conn_pool->acquire();
        try {
            pqxx::work txn(*conn);
            std::cout << graph_name << std::endl;
            pqxx::result rows = txn.exec_params("SELECT * FROM create_graph($1)", graph_name);
            std::cout << "create graph " << graph_name << std::endl;
            txn.commit();
        } catch (const std::exception &e) {
            if ( strstr(e.what(), "already exists") == nullptr ) {
                ERROR(e.what());
                throw std::runtime_error(e.what());
            }
        }
        
        _graph_name = graph_name;
    }

    std::shared_ptr<pqxx::connection> AGEDriver::get_conn() const {
        auto conn = _conn_pool->acquire();
        if (!conn->is_open()) {
            LOG(ERROR, "Database disconnected unexpectedly. ");
            exit(1);
        }
        return conn;
    }
    
    bool AGEDriver::is_directed() const {
        return m_is_directed;
    }

    void AGEDriver::set_worker_thread_num(std::uint64_t new_num) {
        _conn_pool->resize(new_num);
        std::cout<<"set PostgreSQL worker num "<<new_num<<std::endl;
    }

    void AGEDriver::on_edge_writes_finish(){ /* do nothing */ }

    void AGEDriver::thread_exit(){ /* do nothing */ }

    void AGEDriver::on_openmp_workloads_finish() { /* do nothing */ }

    void AGEDriver::configure_distinct_reader_and_writer_threads(std::uint64_t reader_num, std::uint64_t writer_num) { /* do nothing */ }
    
    void AGEDriver::mixed_workload_finish_loading() { /* do nothing */ }

    void AGEDriver::print_and_clear_txn_stats() { /* do nothing */ }

    std::uint64_t AGEDriver::num_vertices() const {
        std::uint64_t result = 0;
        auto conn = get_conn();
        assert(conn->is_open());

        try {
            pqxx::work txn(*conn);
            pqxx::result res = txn.exec_params("SELECT * FROM cypher($1, $$ MATCH (n) RETURN count(n) $$) as (n agtype)", GRAPH_NAME);
            result = res[0]["n"].as<std::uint64_t>();
            txn.commit();
        } catch (const std::exception &e) {
            ERROR(e.what());
            throw std::runtime_error(e.what());
        }

        return result;
    }

    std::uint64_t AGEDriver::num_edges() const {
        std::uint64_t result = 0;
        auto conn = get_conn();
        assert(conn->is_open());

        try {
            pqxx::work txn(*conn);
            pqxx::result res = txn.exec_params("SELECT * FROM cypher($1, $$ MATCH ()-[e]->() RETURN count(e) $$) as (e agtype)", GRAPH_NAME);
            if (res.size() > 0) {
                result = res[0]["e"].as<std::uint64_t>();
            }
        } catch (const std::exception &e) {
            ERROR(e.what());
            throw std::runtime_error(e.what());
        }

        return result;
    }

    std::uint64_t AGEDriver::get_vertex_out_degree(std::uint64_t vertex_id) const {
        std::uint64_t result = 0;
        auto conn = get_conn();

        try {
            pqxx::work txn(*conn);
            pqxx::result res = txn.exec_params("SELECT * FROM cypher($1, $$ MATCH (n {extid: $2})-[e]->() RETURN count(e) $$) as (e agtype)", GRAPH_NAME, vertex_id);
            if (res.size() > 0) {
                result = res[0]["e"].as<std::uint64_t>();
            }
        } catch(const std::exception &e) {
            ERROR(e.what());
            throw std::runtime_error(e.what());
        }
        
        return result;
    }

    std::uint64_t AGEDriver::get_vertex_in_degree(std::uint64_t vertex_id) const {
        std::uint64_t result = 0;
        auto conn = get_conn();

        try {
            pqxx::work txn(*conn);
            pqxx::result res = txn.exec_params("SELECT * FROM cypher($1, $$ MATCH ()-[e]->(n {extid: $2}) RETURN count(e) $$) as (e agtype)", GRAPH_NAME, vertex_id);
            if (res.size() > 0) {
                result = res[0]["e"].as<std::uint64_t>();
            }
        } catch(const std::exception &e) {
            ERROR(e.what());
            throw std::runtime_error(e.what());
        }

        return result;
    }

    vector<std::uint64_t> AGEDriver::get_vertex_out_neighbors(std::uint64_t vertex_id) const {
        std::shared_ptr<pqxx::connection> conn = get_conn();
        assert(conn->is_open());
        std::vector<std::uint64_t> neighbors;

        try {
            pqxx::work txn(*conn);
            pqxx::result res = txn.exec_params(R"sql(SELECT * FROM cypher($1, $$
                    MATCH (n {extid: $2})-[e]->(m) 
                    RETURN m.extid AS o $$) AS (o agtype))sql", GRAPH_NAME, vertex_id);
            txn.commit();
            #pragma omp parallel for
            for (auto const &row : res) {
                neighbors.push_back(row["o"].as<std::uint64_t>());
            }
        } catch (const std::exception &e) {
            ERROR(e.what());
        }

        return neighbors;
    }

    vector<std::uint64_t> AGEDriver::get_vertex_in_neighbors(std::uint64_t vertex_id) const {
        std::shared_ptr<pqxx::connection> conn = get_conn();
        assert(conn->is_open());
        std::vector<std::uint64_t> neighbors;

        try {
            pqxx::work txn(*conn);
            pqxx::result res = txn.exec_params(R"sql(SELECT * FROM cypher($1, $$
                    MATCH (n)-[e]->(m {extid: $2}) 
                    RETURN n.extid AS o $$) AS (o agtype))sql", GRAPH_NAME, vertex_id);
            txn.commit();
            #pragma omp parallel for
            for (auto const &row : res) {
                neighbors.push_back(row["o"].as<std::uint64_t>());
            }
        } catch (const std::exception &e) {
            ERROR(e.what());
        }

        return neighbors;
    }

    std::vector<gfe::graph::WeightedEdge> AGEDriver::get_out_weighted_edges(std::uint64_t vertex_id) const {
        std::shared_ptr<pqxx::connection> conn = get_conn();
        assert(conn->is_open());
        pqxx::work txn(*conn);
        std::vector<gfe::graph::WeightedEdge> neighbors;
        
        try {
            pqxx::result res = txn.exec_params(R"sql(SELECT * FROM cypher($1, $$
                    MATCH (n {extid: $2})-[e]->(m) 
                    RETURN m.extid AS o, e.weight AS w $$) AS (o agtype, w agtype))sql", GRAPH_NAME, vertex_id);
            txn.commit();
            #pragma omp parallel for
            for (auto row : res) {
                double weight = row["w"].as<double>();
                std::uint64_t dest = row["o"].as<std::uint64_t>();
                neighbors.push_back({vertex_id, dest, weight});
            }
        } catch (const std::exception &e) {
            ERROR(e.what());
        }

        return neighbors;
    }

    void AGEDriver::set_timeout(std::uint64_t seconds) {
        m_timeout = chrono::seconds{ seconds };
    }

    void AGEDriver::analytical_workload_end() { /* do nothing */ }

    void AGEDriver::finish_loading() { /* do nothing */ }

    bool AGEDriver::add_vertex(std::uint64_t vertex_id) {
        /* Firstly we check if the same external_id already exists in the graph */
        if (has_vertex(vertex_id)) {
            return false;
        }

        auto conn = get_conn();
        try {
            pqxx::work txn(*conn);
            pqxx::result res = txn.exec_params("SELECT * FROM cypher($1, $$ CREATE (v {extid: $2}) $$) as (v agtype)", GRAPH_NAME, vertex_id);
            txn.commit();

            if (vertex_id > get_max_vertex_extid()) {
                set_max_vertex_extid(vertex_id);
            }
        } catch (const std::exception &e) {
            ERROR(e.what());
            throw std::runtime_error(e.what());
        }
        
        return true;
    }

    // Deleting vertex in age should also delete edges that start or end on said vertex
    bool AGEDriver::remove_vertex(std::uint64_t vertex_id) {
        if (!has_vertex(vertex_id)) {
            return false;
        }

        auto conn = get_conn();
        try {
            pqxx::work txn(*conn);
            pqxx::result res = txn.exec_params("SELECT * FROM cypher($1, $$ MATCH (v {extid: $2}) DETACH DELETE v $$) as (id agtype)", GRAPH_NAME, vertex_id);
            txn.commit();

            if (vertex_id == get_max_vertex_extid()) {
                set_max_vertex_extid(INVALID_EXTID);
            }
        } catch (const std::exception &e) {
            ERROR(e.what());
            throw std::runtime_error(e.what());
        }

        return true;
    }

    bool AGEDriver::has_vertex(std::uint64_t vertex_id) const {
        auto conn = get_conn();
        try {
            pqxx::work txn(*conn);
            pqxx::result res = txn.exec_params(R"sql(SELECT * FROM cypher($1, $$ MATCH (n {extid: $2}) RETURN n $$) as (n agtype))sql", GRAPH_NAME, vertex_id);
            txn.commit();
            return res.size() > 0;
        } catch (const std::exception &e) {
            ERROR(e.what());
            throw std::runtime_error(e.what());
        }

        return false;
    }

    double AGEDriver::get_weight(std::uint64_t source, std::uint64_t destination) const {
        double weight = 0;
        auto conn = get_conn();
        try {
            pqxx::work txn(*conn);
            pqxx::result res = txn.exec_params("SELECT * FROM cypher($1, $$ MATCH (n {extid: $2})-[e: WeightedEdge]->(m {extid: $3}) RETURN e.weight $$) as (weight agtype)", GRAPH_NAME, source, destination);
            txn.commit();
            if (res.size() == 0) {
                weight = std::numeric_limits<double>::quiet_NaN();
            } else if (res.size() > 1) {
                ERROR("More than 1 edge between " << source << " and " << destination);
                throw std::runtime_error("More than 1 edge between " + std::to_string(source) + " and " + std::to_string(destination));
            }

            weight = res[0]["weight"].as<double>();
        } catch (const std::exception &e) {
            ERROR(e.what());
            throw std::runtime_error(e.what());
        }

        return weight;
    }

    bool AGEDriver::add_edge(gfe::graph::WeightedEdge e) {
        vertex_t source_id = e.source();
        vertex_t dest_id = e.destination();
        double weight = e.weight();

        // check whether weighted edge already exists
        if (std::isnan(get_weight(source_id, dest_id))) {
            return false;
        }

        /*
         * To create an edge between two vertices, we first MATCH the two vertices. 
         * Once the nodes are matched, we create an edge between them.
         * 
         *  SELECT * 
         *  FROM cypher('graph_name', $$
         *      MATCH (a {extid: $source_id}), (b {extid: $dest_id})
         *      CREATE (a)-[e:WeightedEdge {weight: $weight}]->(b)
         *      CREATE (b)-[e:WeightedEdge {weight: $weight}]->(a)
         *      RETURN e
         *  $$) as (e agtype);
         * 
         * We don't have to explicitly check if source vertex and destination vertex exist, 
         * since AGE will return the created edge.
         */
        auto conn = get_conn();
        try {
            pqxx::work txn(*conn);
            if (is_directed()) {
                pqxx::result res = txn.exec_params("SELECT * FROM cypher($1, $$ MATCH (a {extid: $2}), (b {extid: $3}) CREATE (a)-[e:WeightedEdge {weight: $4}]->(b) RETURN e $$) as (e agtype)", GRAPH_NAME, source_id, dest_id, weight);
            } else {
                pqxx::result res = txn.exec_params("SELECT * FROM cypher($1, $$ MATCH (a {extid: $2}), (b {extid: $3}) CREATE (a)-[e:WeightedEdge {weight: $4}]->(b) CREATE (b)-[e:WeightedEdge {weight: $4}]->(a) RETURN e $$) as (e agtype)", GRAPH_NAME, source_id, dest_id, weight);
            }
            txn.commit();
        } catch (const std::exception &e) {
            ERROR(e.what());
            throw std::runtime_error(e.what());
        }

        return true;
    }

    bool AGEDriver::add_edge_v2(gfe::graph::WeightedEdge e) {
        vertex_t source = e.source();
        vertex_t destination = e.destination();

        if (!has_vertex(source)) {
            add_vertex(source);
        }
        if (!has_vertex(destination)) {
            add_vertex(destination);
        }
        return add_edge(e);
    }
    
    bool AGEDriver::add_edge_v3(gfe::graph::WeightedEdge edge){
        vertex_t source = edge.source();
        vertex_t destination = edge.destination();
        double weight = edge.weight();
        double prev_weight = get_weight(source, destination);

        if (!std::isnan(prev_weight) && !DBL_EQ(prev_weight, weight)) {
            auto conn = get_conn();
            try {
                pqxx::work txn(*conn);
                if (is_directed()) {
                    pqxx::result res = txn.exec_params("SELECT * FROM cypher($1, $$ MATCH (n {extid: $2})-[e: WeightedEdge {weight: $3}]->(m {extid: $4}) SET e.weight=$5 $$) as (e agtype)", GRAPH_NAME, source, prev_weight, destination, weight);
                } else {
                    pqxx::result res = txn.exec_params("SELECT * FROM cypher($1, $$ MATCH (a)-[r1]->(b), (b)-[r2]->(a) WHERE a.extid = $2 AND b.extid = $3 AND r1.weight = $4 AND r2.weight = $4 SET r1.weight=$5, r2.weight=$5 $$) as (e agtype)", GRAPH_NAME, source, destination, prev_weight, weight);
                }
                txn.commit();
            } catch (const std::exception &e) {
                ERROR(e.what());
                throw std::runtime_error(e.what());
            }
        } else {
            // create edge
            add_edge(edge);
        }

        // If the edge not exists, insert a new edge
        return add_edge_v2(edge);
    }

    /*
     * lookup the edge, and then insert/update with (old_weight + new_weight)
     */
    bool AGEDriver::update_edge_v1(gfe::graph::WeightedEdge edge) {
        vertex_t source = edge.source();
        vertex_t destination = edge.destination();
        double weight = edge.weight();

        double prev_weight = get_weight(source, destination);
        if (!std::isnan(prev_weight)) {
            // If the edge already exists, its weight is updated
            weight += prev_weight;
            add_edge_v3({source, destination, weight});
        } else {
            // create edge
            add_edge_v3(edge);
        }

        return true;
    }

    bool AGEDriver::remove_edge(gfe::graph::Edge e){
       vertex_t source = e.source();
       vertex_t destination = e.destination();

        auto conn = get_conn();
        try {
            pqxx::work txn(*conn);
            if (is_directed()) {
                pqxx::result res = txn.exec_params("SELECT * FROM cypher($1, $$ MATCH (n {extid: $2})-[r]->(m {extid: $3}) DELETE r $$) as (e agtype)", GRAPH_NAME, source, destination);
            } else {
                pqxx::result res = txn.exec_params("SELECT * FROM cypher($1, $$ MATCH (a)-[r1]->(b), (b)-[r2]->(a) WHERE a.extid = $2 AND b.extid = $3 DELETE r1, r2 $$) as (e agtype)", GRAPH_NAME, source, destination);
            }
            txn.commit();
        } catch (const std::exception &e) {
            ERROR(e.what());
            throw std::runtime_error(e.what());
        }

        return true;
    }


    /*****************************************************************************
    *                                                                           *
    *  Dump                                                                     *
    *                                                                           *
    *****************************************************************************/
    void AGEDriver::dump_ostream(std::ostream& out) const {
        out << "[AGE] num vertices: " << m_num_vertices << ", num edges: " << m_num_edges << ", "
                                                                                                   "directed graph: " << boolalpha << is_directed() << ", read only txn for graphalytics: " << m_read_only << endl;
        const std::uint64_t max_vertex_id = get_max_vertex_extid();
        for(std::uint64_t external_id = 1; external_id <= max_vertex_id; external_id++){
            if(!has_vertex(external_id)) continue; // the vertex has been deleted
            out << "[external_id: " << external_id << "]";
            { // outgoing edges
                out << " outgoing edges: ";
                auto neighbors = get_out_weighted_edges(external_id);
                bool first = true;
                for (auto edge : neighbors){
                    if(first){ first = false; } else { out << ", "; }
                    std::uint64_t dst_id = edge.destination();
                    double weight = edge.weight();
                    out << "<" << " [extid: " << dst_id << "], " << weight << ">";
                }
            }
            out << endl;
        }
    }

    /*****************************************************************************
     *                                                                           *
     *  BFS                                                                      *
     *                                                                           *
     *****************************************************************************/
    // Implementation based on the reference BFS for the GAP Benchmark Suite
    // https://github.com/sbeamer/gapbs
    // The reference implementation has been written by Scott Beamer
    //
    // Copyright (c) 2015, The Regents of the University of California (Regents)
    // All Rights Reserved.
    //
    // Redistribution and use in source and binary forms, with or without
    // modification, are permitted provided that the following conditions are met:
    // 1. Redistributions of source code must retain the above copyright
    //    notice, this list of conditions and the following disclaimer.
    // 2. Redistributions in binary form must reproduce the above copyright
    //    notice, this list of conditions and the following disclaimer in the
    //    documentation and/or other materials provided with the distribution.
    // 3. Neither the name of the Regents nor the
    //    names of its contributors may be used to endorse or promote products
    //    derived from this software without specific prior written permission.
    //
    // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
    // ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    // WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    // DISCLAIMED. IN NO EVENT SHALL REGENTS BE LIABLE FOR ANY
    // DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    // (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    // LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    // ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    // SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    /*

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
    //#define DEBUG_BFS
#if defined(DEBUG_BFS)
#define COUT_DEBUG_BFS(msg) COUT_DEBUG(msg)
#else
#define COUT_DEBUG_BFS(msg)
#endif

    void AGEDriver::bfs(std::uint64_t external_source_id, const char* dump2file) {
        //std::cout<<"from source "<<external_source_id<<std::endl;
        if(m_is_directed) { ERROR("This implementation of the BFS does not support directed graphs"); }

        // Init
        utility::TimeoutService timeout { m_timeout };
        Timer timer; timer.start();
        std::uint64_t max_vertex_id = get_max_vertex_extid();
        std::unique_ptr<int64_t[]> ptr_distances;
        //std::uint64_t num_edges = m_num_edges;
        //COUT_DEBUG_BFS("root: " << root << " [external vertex: " << external_source_id << "]");

        // Run the BFS algorithm
        auto conn = get_conn();
        try {
            pqxx::work txn(*conn);
            pqxx::result res = txn.exec_params(R"sql(SELECT * FROM cypher($1, $$
                    WITH RECURSIVE bfs AS (
                        SELECT id, extid, 1 AS level, ARRAY[id] AS visited
                        FROM cypher($1, $$
                            MATCH (n {extid: $2})
                            RETURN id(n) AS id, n.extid AS extid
                        $$) AS (id agtype, extid agtype)
                        
                        UNION ALL
                        
                        SELECT e.end_id AS id, e.end_extid AS extid, p.level + 1 AS level, p.visited || e.end_id AS visited
                        FROM bfs AS p
                        JOIN cypher($3, $$
                            MATCH (n)-[r]->(m)
                            RETURN id(n) AS start_id, id(m) AS end_id, m.extid AS end_extid
                        $$) AS e(start_id agtype, end_id agtype, end_extid agtype) ON e.start_id = p.id
                        WHERE e.end_id<>ALL(p.visited)
                    )
                    SELECT extid, level FROM bfs
                $$) as (extid agtype, level agtype))sql", GRAPH_NAME, external_source_id, GRAPH_NAME);
            txn.commit();

            #pragma omp parallel
            {
                #pragma omp for
                for (std::uint64_t i = 0; i < max_vertex_id; i++) {
                    ptr_distances[i] = std::numeric_limits<std::uint64_t>::max();
                }

                #pragma omp for
                for (size_t i = 0; i < res.size(); i++) {
                    std::uint64_t vertex_id = res[i]["extid"].as<std::uint64_t>();
                    std::uint64_t distance = res[i]["level"].as<std::uint64_t>();
                    ptr_distances[vertex_id] = distance;
                }
            }
        } catch (const std::exception &e) {
            ERROR(e.what());
            exit(1);
        }

        //unique_ptr<int64_t[]> ptr_result = static_do_bfs(transaction, num_vertices, num_edges, max_vertex_id, root, timeout);
        //cout << "BFS took " << t << endl;
        if(timeout.is_timeout()){
            RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);
        }

        std::vector<std::pair<std::uint64_t, std::int64_t>> results;
        for (std::uint64_t i = 0; i < max_vertex_id; i++) {
            results.push_back({i, ptr_distances[i]});
        }

        //cout << "Translation took " << t << endl;
        //std::cout<<dump2file<<std::endl;
        if(dump2file != nullptr) // store the results in the given file
            save_results<int64_t, false>(results, dump2file);
        //std::cout<<"bfs over"<<std::endl;
    }

    /*****************************************************************************
    *                                                                           *
    *  Graphalytics Helpers                                                     *
    *                                                                           *
    *****************************************************************************/

    std::uint64_t AGEDriver::get_max_vertex_extid() const {
        if (_max_vertex_extid.load() == INVALID_EXTID) {
            // find the maximum vertex id
            auto conn = get_conn();
            try {
                pqxx::work txn(*conn);
                pqxx::result res = txn.exec_params(R"sql(SELECT * FROM cypher($1, $$
                        MATCH (n) 
                        RETURN max(n.extid) 
                    $$) AS (max_extid agtype))sql", GRAPH_NAME);
                txn.commit();
                if (res.size() > 0) {
                    _max_vertex_extid.store(res[0]["max_extid"].as<std::uint64_t>());
                }
            } catch (const std::exception &e) {
                ERROR(e.what());
                throw std::runtime_error(e.what());
            }
        }

        return _max_vertex_extid.load();
    }

    void AGEDriver::set_max_vertex_extid(std::uint64_t max_vertex_extid) {
        _max_vertex_extid.store(max_vertex_extid);
    }
    
    template <typename T, bool negative_scores>
    void AGEDriver::save_results(const vector<pair<std::uint64_t, T>>& result, const char* dump2file) {
        assert(dump2file != nullptr);
        COUT_DEBUG("save the results to: " << dump2file);

        fstream handle(dump2file, ios_base::out);
        if (!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

        for (const auto &p : result) {
            if(p.first == numeric_limits<std::uint64_t>::max()) continue; // invalid node

            handle << p.first << " ";

            if(!negative_scores && p.second < 0){
                handle << numeric_limits<T>::max();
            } else {
                handle << p.second;
            }

            handle << "\n";
        }

        handle.close();
    }

    /*****************************************************************************
 *                                                                           *
 *  PageRank                                                                 *
 *                                                                           *
 *****************************************************************************/
    //#define DEBUG_PAGERANK
#if defined(DEBUG_PAGERANK)
#define COUT_DEBUG_PAGERANK(msg) COUT_DEBUG(msg)
#else
#define COUT_DEBUG_PAGERANK(msg)
#endif

// Implementation based on the reference PageRank for the GAP Benchmark Suite
// https://github.com/sbeamer/gapbs
// The reference implementation has been written by Scott Beamer
//
// Copyright (c) 2015, The Regents of the University of California (Regents)
// All Rights Reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
// 3. Neither the name of the Regents nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL REGENTS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

/*
GAP Benchmark Suite
Kernel: PageRank (PR)
Author: Scott Beamer

Will return pagerank scores for all vertices once total change < epsilon

This PR implementation uses the traditional iterative approach. This is done
to ease comparisons to other implementations (often use same algorithm), but
it is not necessarily the fastest way to implement it. It does perform the
updates in the pull direction to remove the need for atomics.
*/
    void AGEDriver::pagerank(std::uint64_t num_iterations, double damping_factor, const char* dump2file) {
        if(m_is_directed) { ERROR("This implementation of PageRank does not support directed graphs"); }

        // Init
        utility::TimeoutService timeout { m_timeout };
        Timer timer; timer.start();
        //gt::SharedROTransaction transaction = GTX->begin_shared_read_only_transaction();
        std::uint64_t num_vertices = this->num_vertices();
        std::uint64_t max_vertex_id = get_max_vertex_extid();
        
        // Run the PageRank algorithm
        const double init_score = 1.0 / num_vertices;
        const double base_score = (1.0 - damping_factor) / num_vertices;

        unique_ptr<double[]> ptr_scores{ new double[max_vertex_id]() }; // avoid memory leaks
        unique_ptr<std::uint64_t[]> ptr_degrees{ new std::uint64_t[max_vertex_id]() }; // avoid memory leaks
        double* scores = ptr_scores.get();
        std::uint64_t* __restrict degrees = ptr_degrees.get();

        for (std::uint64_t v = 1; v <= max_vertex_id; v++) {
            scores[v - 1] = init_score;

            // compute the outdegree of the vertex
            if (!has_vertex(v)) { // check the vertex exists
                degrees[v - 1] = get_vertex_out_degree(v);
            } else {
                degrees[v - 1] = numeric_limits<std::uint64_t>::max();
            }
        }

        gapbs::pvector<double> outgoing_contrib(max_vertex_id, 0.0);

        // pagerank iterations
        for(std::uint64_t iteration = 0; iteration < num_iterations && !timeout.is_timeout(); iteration++){
            double dangling_sum = 0.0;

            // for each node, precompute its contribution to all of its outgoing neighbours and, if it's a sink,
            // add its rank to the `dangling sum' (to be added to all nodes).
#pragma omp parallel for reduction(+:dangling_sum)
            for(std::uint64_t v = 1; v <= max_vertex_id; v++){
                std::uint64_t out_degree = degrees[v-1];
                if(out_degree == numeric_limits<std::uint64_t>::max()){
                    continue; // the vertex does not exist
                } else if (out_degree == 0){ // this is a sink
                    dangling_sum += scores[v-1];
                } else {
                    outgoing_contrib[v-1] = scores[v-1] / out_degree;
                }
            }

            dangling_sum /= num_vertices;

            // compute the new score for each node in the graph
            //fixme: currently our txn does not support being executed by mutiple threads. It requires fetching thread_id for each operation execution from the table.
            for (std::uint64_t v = 1; v <= max_vertex_id; v++) {
                if (degrees[v-1] == numeric_limits<std::uint64_t>::max()) { continue; } // the vertex does not exist

                double incoming_total = 0;
                //auto iterator = transaction.simple_get_edges(v, /* label ? */1,thread_id); // fixme: incoming edges for directed graphs
                auto out_neighbors = get_vertex_out_neighbors(v);
                for (auto neighbor : out_neighbors) {
                    incoming_total += outgoing_contrib[neighbor-1];
                }
                // update the score
                scores[v-1] = base_score + damping_factor * (incoming_total + dangling_sum);
            }
        }

        std::vector<std::pair<std::uint64_t, double>> results(max_vertex_id);
        #pragma omp parallel for
        for (std::uint64_t i = 0; i < max_vertex_id; i++) {
            results.push_back({i, scores[i]});
        }
        #pragma omp barrier

        if(dump2file != nullptr)
            save_results(results,dump2file);

        if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);  }
    }

/*****************************************************************************
 *                                                                           *
 *  WCC                                                                      *
 *                                                                           *
 *****************************************************************************/
// Implementation based on the reference WCC for the GAP Benchmark Suite
// https://github.com/sbeamer/gapbs
// The reference implementation has been written by Scott Beamer
//
// Copyright (c) 2015, The Regents of the University of California (Regents)
// All Rights Reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
// 3. Neither the name of the Regents nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL REGENTS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#define DEBUG_WCC
#if defined(DEBUG_WCC)
#define COUT_DEBUG_WCC(msg) COUT_DEBUG(msg)
#else
#define COUT_DEBUG_WCC(msg)
#endif
/*
GAP Benchmark Suite
Kernel: Connected Components (CC)
Author: Scott Beamer

Will return comp array labelling each vertex with a connected component ID

This CC implementation makes use of the Shiloach-Vishkin [2] algorithm with
implementation optimizations from Bader et al. [1]. Michael Sutton contributed
a fix for directed graphs using the min-max swap from [3], and it also produces
more consistent performance for undirected graphs.

[1] David A Bader, Guojing Cong, and John Feo. "On the architectural
    requirements for efficient execution of graph algorithms." International
    Conference on Parallel Processing, Jul 2005.

[2] Yossi Shiloach and Uzi Vishkin. "An o(logn) parallel connectivity algorithm"
    Journal of Algorithms, 3(1):57–67, 1982.

[3] Kishore Kothapalli, Jyothish Soman, and P. J. Narayanan. "Fast GPU
    algorithms for graph connectivity." Workshop on Large Scale Parallel
    Processing, 2010.
*/

// The hooking condition (comp_u < comp_v) may not coincide with the edge's
// direction, so we use a min-max swap such that lower component IDs propagate
// independent of the edge's direction.
/*
    Should we use tarjan algorithm for WCC?
*/
    void AGEDriver::wcc(const char* dump2file) {
        utility::TimeoutService timeout { m_timeout };
        Timer timer; timer.start();
        std::uint64_t max_vertex_id = get_max_vertex_extid();

        // run wcc
        // init
        COUT_DEBUG_WCC("max_vertex_id: " << max_vertex_id);
        unique_ptr<std::uint64_t[]> ptr_components { new std::uint64_t[max_vertex_id] };
        std::uint64_t* comp = ptr_components.get();

#pragma omp parallel for
        for (std::uint64_t n = 1; n <= max_vertex_id; n++){
            if(has_vertex(n)){ // the vertex does not exist
                COUT_DEBUG_WCC("Vertex #" << n << " does not exist");
                comp[n] = numeric_limits<std::uint64_t>::max();
            } else {
                comp[n] = n;
            }
        }

        bool change = true;
        while (change && !timeout.is_timeout()) {
            change = false;

#pragma omp parallel for schedule(dynamic, 64)
            for (std::uint64_t u = 1; u <= max_vertex_id; u++){
                if(comp[u] == numeric_limits<std::uint64_t>::max()) continue; // the vertex does not exist

                auto neighbors = get_vertex_out_neighbors(u);
                for (auto v : neighbors) {
                    std::uint64_t comp_u = comp[u];
                    std::uint64_t comp_v = comp[v];
                    if (comp_u != comp_v) {
                        // Hooking condition so lower component ID wins independent of direction
                        std::uint64_t high_comp = std::max(comp_u, comp_v);
                        std::uint64_t low_comp = std::min(comp_u, comp_v);
                        if (high_comp == comp[high_comp]) {
                            change = true;
                            COUT_DEBUG_WCC("comp[" << high_comp << "] = " << low_comp);
                            comp[high_comp] = low_comp;
                        }
                    }

                }
            }

#pragma omp parallel for schedule(dynamic, 64)
            for (std::uint64_t n = 1; n <= max_vertex_id; n++){
                if(comp[n] == numeric_limits<std::uint64_t>::max()) continue; // the vertex does not exist

                while (comp[n] != comp[comp[n]]) {
                    comp[n] = comp[comp[n]];
                }
            }


            COUT_DEBUG_WCC("change: " << change);
        }

        if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

        std::vector<std::pair<std::uint64_t, std::uint64_t>> results(max_vertex_id);
        #pragma omp parallel for
        for (std::uint64_t i = 0; i < max_vertex_id; i++) {
            results[i] = {i, comp[i]};
        }
        #pragma omp barrier

        // store the results in the given file
        if(dump2file != nullptr)
            save_results(results, dump2file);
        
        // do_weighted_scan(transaction,max_vertex_id);
    }

/*****************************************************************************
 *                                                                           *
 *  CDLP                                                                     *
 *                                                                           *
 *****************************************************************************/
    void AGEDriver::generate_two_hops_neighbor_candidates(std::vector<std::uint64_t>&vertices){
        vertices.resize(two_hop_neighbor_size);
        std::uint64_t max_vertex_id = get_max_vertex_extid();
        std::unordered_set<std::uint64_t> unique_vertices;
        for(std::uint64_t i = 0; i < two_hop_neighbor_size; i++){
            std::uint64_t vid = rand() % max_vertex_id + 1;
            auto result = unique_vertices.emplace(vid);
            if(result.second){
                vertices[i]=vid;
            }else{
                i--;
            }
        }
    }

    void AGEDriver::cdlp(std::uint64_t max_iterations, const char* dump2file) {
     /*   if(m_is_directed) { ERROR("This implementation of the CDLP does not support directed graphs"); } */

        utility::TimeoutService timeout { m_timeout };
        Timer timer; timer.start();
        std::uint64_t max_vertex_id = get_max_vertex_extid();

        // Run the CDLP algorithm
        unique_ptr<std::uint64_t[]> ptr_labels0 { new std::uint64_t[max_vertex_id] };
        unique_ptr<std::uint64_t[]> ptr_labels1 { new std::uint64_t[max_vertex_id] };
        std::uint64_t* labels0 = ptr_labels0.get(); // current labels
        std::uint64_t* labels1 = ptr_labels1.get(); // labels for the next iteration

        // initialisation
#pragma omp parallel for
        for(std::uint64_t v = 1; v <= max_vertex_id; v++){
            if(has_vertex(v)){ // the vertex does not exist
                labels0[v] = labels1[v] = numeric_limits<std::uint64_t>::max();
            } else {
                labels0[v] = v;
            }
        }

        // algorithm pass
        bool change = true;
        std::uint64_t current_iteration = 0;
        while(current_iteration < max_iterations && change && !timeout.is_timeout()){
            change = false; // reset the flag

#pragma omp parallel for schedule(dynamic, 64) shared(change)
            for(std::uint64_t v = 1; v <= max_vertex_id; v++){
                if(labels0[v] == numeric_limits<std::uint64_t>::max()) continue; // the vertex does not exist

                unordered_map<std::uint64_t, std::uint64_t> histogram;

                // compute the histogram from both the outgoing & incoming edges. The aim is to find the number of each label
                // is shared among the neighbours of node_id
                vector<std::uint64_t> out_neighbors = get_vertex_out_neighbors(v); // out edges
                for (auto u : out_neighbors) {
                    histogram[labels0[u]]++;
                }
                // get the max label
                std::uint64_t label_max = numeric_limits<int64_t>::max();
                std::uint64_t count_max = 0;
                for(const auto pair : histogram){
                    if(pair.second > count_max || (pair.second == count_max && pair.first < label_max)){
                        label_max = pair.first;
                        count_max = pair.second;
                    }
                }

                labels1[v] = label_max;
                change |= (labels0[v] != labels1[v]);
            }

            std::swap(labels0, labels1); // next iteration
            current_iteration++;
        }

        if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

        // collect results
        std::unique_ptr<std::uint64_t[]> ptr_labels = labels0 == ptr_labels0.get() ? std::move(ptr_labels0) : std::move(ptr_labels1);
        std::vector<std::pair<std::uint64_t, std::uint64_t>> results(max_vertex_id);
        #pragma omp parallel for
        for (std::uint64_t i = 0; i < max_vertex_id; i++) {
            results[i] = {i, ptr_labels[i]};
        }
        #pragma omp barrier

        // Store the results in the given file
        if(dump2file != nullptr)
            save_results(results, dump2file);
    }

    void AGEDriver::one_hop_neighbors(std::vector<std::uint64_t>&vertices){
        std::unordered_map<std::uint64_t, std::vector<std::uint64_t>> results;

         #pragma omp parallel for
         for (auto source : vertices) {
             results[source] = get_vertex_out_neighbors(source);
         }
         //omp for finished
    }

    void AGEDriver::two_hop_neighbors(std::vector<std::uint64_t>&vertices){
        std::unordered_map<std::uint64_t, std::vector<std::uint64_t>> results;

        #pragma omp parallel for
        for (auto source : vertices) {
            try {
                auto conn = get_conn();
                std::vector<std::uint64_t> neighbors;
                // get the two-hop neighbors
                pqxx::work txn(*conn);
                pqxx::result res = txn.exec_params(R"sql(
                SELECT * FROM cypher($1, $$
                    MATCH (n {extid: $2})-[r1]->(m)-[r2]->(o) 
                    RETURN o.extid AS o $$) AS (o agtype))sql", GRAPH_NAME, source);
                txn.commit();

                for (auto row : res) {
                    neighbors.push_back(row["o"].as<std::uint64_t>());
                }
                #pragma omp critical(two_hop_neighbors)
                {
                    results[source] = neighbors;
                }
            } catch (const std::exception &e) {
                ERROR(e.what());
            }
        }
    }

/*****************************************************************************
 *                                                                           *
 *  LCC                                                                      *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG_LCC
#if defined(DEBUG_LCC)
#define COUT_DEBUG_LCC(msg) COUT_DEBUG(msg)
#else
#define COUT_DEBUG_LCC(msg)
#endif
// loosely based on the impl~ made for GraphOne
    void AGEDriver::lcc(const char* dump2file) {
        if(m_is_directed) { ERROR("Implementation of LCC supports only undirected graphs"); }

        // Init
        utility::TimeoutService timeout { m_timeout };
        Timer timer; timer.start();
        std::uint64_t max_vertex_id = get_max_vertex_extid();

        // Run the LCC algorithm
        unique_ptr<double[]> ptr_lcc { new double[max_vertex_id] };
        double* lcc = ptr_lcc.get();
        unique_ptr<uint32_t[]> ptr_degrees_out { new uint32_t[max_vertex_id] };
        uint32_t* __restrict degrees_out = ptr_degrees_out.get();

        // precompute the degrees of the vertices
#pragma omp parallel for schedule(dynamic, 4096)
        for(std::uint64_t v = 1; v <= max_vertex_id; v++){
            if(!has_vertex(v)){
                lcc[v] = numeric_limits<double>::signaling_NaN();
            } else {
                degrees_out[v] = get_vertex_out_degree(v);
            }
        }

#pragma omp parallel for schedule(dynamic, 64)
        for(std::uint64_t v = 1; v <= max_vertex_id; v++){
            if(degrees_out[v] == numeric_limits<uint32_t>::max()) continue; // the vertex does not exist

            COUT_DEBUG_LCC("> Node " << v);
            if(timeout.is_timeout()) continue; // exhausted the budget of available time
            lcc[v] = 0.0;
            std::uint64_t num_triangles = 0; // number of triangles found so far for the node v

            // Cfr. Spec v.0.9.0 pp. 15: "If the number of neighbors of a vertex is less than two, its coefficient is defined as zero"
            std::uint64_t v_degree_out = degrees_out[v];
            if(v_degree_out < 2) continue;

            // Build the list of neighbours of v
            vector<std::uint64_t> edges = get_vertex_out_neighbors(v);
            unordered_set<std::uint64_t> neighbours(edges.begin(), edges.end());

            // again, visit all neighbours of v
            // for directed graphs, edges1 contains the intersection of both the incoming and the outgoing edges
            for (auto u : edges) {
                COUT_DEBUG_LCC("[" << i << "/" << edges.size() << "] neighbour: " << u);
                assert(neighbours.count(u) == 1 && "The set `neighbours' should contain all neighbours of v");

                // For the Graphalytics spec v 0.9.0, only consider the outgoing edges for the neighbours u
                auto neighbor_of_u = get_vertex_out_neighbors(u);
                for (auto w : neighbor_of_u) {
                    COUT_DEBUG_LCC("---> [" << j << "/" << /* degree */ (u_out_interval.second - u_out_interval.first) << "] neighbour: " << w);
                    // check whether it's also a neighbour of v
                    if(neighbours.count(w) == 1){
                        COUT_DEBUG_LCC("Triangle found " << v << " - " << u << " - " << w);
                        num_triangles++;
                    }
                }
            }

            // register the final score
            std::uint64_t max_num_edges = v_degree_out * (v_degree_out -1);
            lcc[v] = static_cast<double>(num_triangles) / max_num_edges;
            COUT_DEBUG_LCC("Score computed: " << (num_triangles) << "/" << max_num_edges << " = " << lcc[v]);
        }

        if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);  }

        std::vector<std::pair<std::uint64_t, double>> results(max_vertex_id);
        #pragma omp parallel for
        for (std::uint64_t i = 0; i < max_vertex_id; i++) {
            results.push_back({i, lcc[i]});
        }
        #pragma omp barrier

        // Store the results in the given file
        if(dump2file != nullptr)
            save_results(results, dump2file);
    }

/*****************************************************************************
 *                                                                           *
 *  SSSP                                                                     *
 *                                                                           *
 *****************************************************************************/
// Implementation based on the reference SSSP for the GAP Benchmark Suite
// https://github.com/sbeamer/gapbs
// The reference implementation has been written by Scott Beamer
//
// Copyright (c) 2015, The Regents of the University of California (Regents)
// All Rights Reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
// 3. Neither the name of the Regents nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL REGENTS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    using NodeID = std::uint64_t;
    using WeightT = double;
    static const size_t kMaxBin = numeric_limits<size_t>::max()/2;
    // static
    // gapbs::pvector<WeightT> do_sssp(gfe::library::AGEDriver* driver, NodeID num_edges, NodeID max_vertex_id, NodeID source, double delta, utility::TimeoutService& timer) {


    //     return dist;
    // }

    void AGEDriver::sssp(std::uint64_t source_vertex_id, const char* dump2file) {
        std::cout<<"sssp running"<<std::endl;

       utility::TimeoutService timeout { m_timeout };
        Timer timer; timer.start();
        NodeID &source = source_vertex_id;
        std::uint64_t num_edges = this->num_edges();
        std::uint64_t max_vertex_id = get_max_vertex_extid();

        // Run the SSSP algorithm
        // Init
        double delta = 2.0; // same value used in the GAPBS, at least for most graphs
        gapbs::pvector<WeightT> dist(max_vertex_id, numeric_limits<WeightT>::infinity());
        dist[source-1] = 0;
        gapbs::pvector<NodeID> frontier(num_edges);
        // two element arrays for double buffering curr=iter&1, next=(iter+1)&1
        size_t shared_indexes[2] = {0, kMaxBin};
        size_t frontier_tails[2] = {1, 0};
        frontier[0] = source;
#pragma omp parallel
        {
            vector<vector<NodeID> > local_bins(0);
            size_t iter = 0;

            while (shared_indexes[iter&1] != kMaxBin) {
                size_t &curr_bin_index = shared_indexes[iter&1];
                size_t &next_bin_index = shared_indexes[(iter+1)&1];
                size_t &curr_frontier_tail = frontier_tails[iter&1];
                size_t &next_frontier_tail = frontier_tails[(iter+1)&1];
#pragma omp for nowait schedule(dynamic, 64)
                for (size_t i=0; i < curr_frontier_tail; i++) {
                    NodeID u = frontier[i];
                    if (dist[u-1] >= delta * static_cast<WeightT>(curr_bin_index)) {
                        std::vector<gfe::graph::WeightedEdge> edges = get_out_weighted_edges(u);
                        for (gfe::graph::WeightedEdge &edge : edges) {
                            std::uint64_t v = edge.destination();
                            double w = edge.weight();
                            WeightT old_dist = dist[v-1];
                            WeightT new_dist = dist[u-1] + w;
                            if (new_dist < old_dist) {
                                bool changed_dist = true;
                                while (!gapbs::compare_and_swap(dist[v-1], old_dist, new_dist)) {
                                    old_dist = dist[v-1];
                                    if (old_dist <= new_dist) {
                                        changed_dist = false;
                                        break;
                                    }
                                }
                                if (changed_dist) {
                                    size_t dest_bin = new_dist/delta;
                                    if (dest_bin >= local_bins.size()) {
                                        local_bins.resize(dest_bin+1);
                                    }
                                    local_bins[dest_bin].push_back(v);
                                }
                            }

                        }
                    }
                }
                for (size_t i=curr_bin_index; i < local_bins.size(); i++) {
                    if (!local_bins[i].empty()) {
#pragma omp critical
                        next_bin_index = min(next_bin_index, i);
                        break;
                    }
                }

#pragma omp barrier
#pragma omp single nowait
                {
                    curr_bin_index = kMaxBin;
                    curr_frontier_tail = 0;
                }

                if (next_bin_index < local_bins.size()) {
                    size_t copy_start = gapbs::fetch_and_add(next_frontier_tail, local_bins[next_bin_index].size());
                    copy(local_bins[next_bin_index].begin(), local_bins[next_bin_index].end(), frontier.data() + copy_start);
                    local_bins[next_bin_index].resize(0);
                }

                iter++;
#pragma omp barrier
            }

#if defined(DEBUG)
            #pragma omp single
        COUT_DEBUG("took " << iter << " iterations");
#endif
        }
        // auto distances = do_static_sssp(num_edges, max_vertex_id, root, delta, timeout);
        if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);  }

        std::vector<std::pair<std::uint64_t, double>> results(max_vertex_id);
        #pragma omp parallel for
        for (std::uint64_t i = 0; i < max_vertex_id; i++) {
            results.push_back({i, dist[i]});
        }
        #pragma omp barrier

        // Store the results in the given file
        if(dump2file != nullptr)
            save_results(results, dump2file);
        //std::cout<<"sssp ran"<<std::endl;
    }
}//namespace gfe::library
