//
// Created by Ziming Zhang on 8/13/24.
//

#pragma once
#include <atomic>
#include <chrono>
#include <vector>
#include <string>
#include <queue>
#include <memory>
#include <mutex>
#include <condition_variable>
#include <limits>
#include <float.h>
#include <cstdint>
//#include "library/interface.hpp"
#include "../interface.hpp"
#include "../../graph/edge.hpp"
#include <tbb/enumerable_thread_specific.h>//to count the time
#include <pqxx/pqxx>

#define GTX_SET_THREAD_NUM true

namespace gfe::library {
#define COUT_DEBUG_FORCE(msg) { std::scoped_lock<std::mutex> lock{::gfe::_log_mutex}; std::cout << "[TeseoDriver::" << __FUNCTION__ << "] " << msg << std::endl; }
#if defined(DEBUG)
#define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
#define COUT_DEBUG(msg)
#endif
    
    class ConnectionPool;

    class AGEDriver: public virtual UpdateInterface, public virtual GraphalyticsInterface {
    private: 

        AGEDriver(const AGEDriver&)= delete;
        AGEDriver& operator=(const AGEDriver&)= delete;

    protected:
        const bool m_is_directed; // whether the underlying graph is directed or undirected
        const bool m_read_only; // whether to used read only transactions for graphalytics
        mutable std::atomic<uint64_t> _max_vertex_extid = INVALID_EXTID;
        std::atomic<std::uint64_t> m_num_vertices {0}; // keep track of the total number of vertices
        std::atomic<std::uint64_t> m_num_edges {0}; // keep track of the total number fo edges
        std::chrono::seconds m_timeout {0}; // the budget to complete each of the algorithms in the Graphalytics suite

        std::uint64_t get_max_vertex_extid() const;

        void set_max_vertex_extid(std::uint64_t max_vertex_extid);

        // Helper, save the content of the vector to the given output file
        template <typename T, bool negative_scores = true>
        void save_results(const std::vector<std::pair<std::uint64_t, T>>& result, const char* dump2file);
        
        /*
         * connect to a postgresql backend
         */
        void connect_to_database(size_t pool_size, const std::string& host, const std::string& port, const std::string& dbname, const std::string& user, const std::string& password);
        void open_graph(const std::string& graph_name);

    public:
        std::uint64_t get_vertex_out_degree(std::uint64_t vertex_id) const;

        std::uint64_t get_vertex_in_degree(std::uint64_t vertex_id) const;

        std::vector<std::uint64_t> get_vertex_out_neighbors(std::uint64_t vertex_id) const;

        std::vector<std::uint64_t> get_vertex_in_neighbors(std::uint64_t vertex_id) const;

        /* return the weighted edges of the vertex */
        std::vector<gfe::graph::WeightedEdge> get_out_weighted_edges(std::uint64_t vertex_id) const;

    public:
        AGEDriver(bool is_directed, bool read_only = true);

        virtual void set_worker_thread_num(std::uint64_t new_num);
        virtual void on_edge_writes_finish();
        virtual void finish_loading();
        virtual void thread_exit();
        virtual void on_openmp_workloads_finish();
        virtual void configure_distinct_reader_and_writer_threads(std::uint64_t reader_num, std::uint64_t writer_num);
        //debug information
        virtual void print_and_clear_txn_stats();
        virtual void analytical_workload_end();

        virtual void mixed_workload_finish_loading();

        virtual ~AGEDriver() noexcept ;

        virtual std::uint64_t num_edges() const;

        virtual std::uint64_t num_vertices() const;

        /**
         * Returns true if the given vertex is present, false otherwise
         */
        virtual bool has_vertex(std::uint64_t vertex_id) const;

        /**
         * Returns the weight of the given edge is the edge is present, or NaN otherwise
         */
        virtual double get_weight(std::uint64_t source, std::uint64_t destination) const;

        /**
         * Check whether the graph is directed
         */
        virtual bool is_directed() const;

        /**
         * Impose a timeout on each graph computation. A computation that does not terminate by the given seconds will raise a TimeoutError.
         */
        virtual void set_timeout(std::uint64_t seconds);

        /**
         * Add the given vertex to the graph
         * @return true if the vertex has been inserted, false otherwise (that is, the vertex already exists)
         */
        virtual bool add_vertex(std::uint64_t vertex_id);
        /**
     * Remove the mapping for a given vertex. The actual internal vertex is not removed from the adjacency list.
     * @param vertex_id the vertex to remove
     * @return true if a mapping for that vertex existed, false otherwise
     */
        virtual bool remove_vertex(std::uint64_t vertex_id);

        /**
         * Add the given edge in the graph if it doesn't exist
         * @return true if the edge has been inserted, false if this edge already exists or one of the referred
         *         vertices does not exist.
         */
        virtual bool add_edge(gfe::graph::WeightedEdge e);

        /**
        * Add the given edge in the graph. Implicitly create the referred vertices if they do not already exist.
        * If the edge already exists, its weight is updated.
        * @return always true.
        */
        virtual bool add_edge_v2(gfe::graph::WeightedEdge e);
        /**
         * same as v2 for age
         * @param e
         * @return
         */
        virtual bool add_edge_v3(gfe::graph::WeightedEdge e);

        virtual bool update_edge_v1(gfe::graph::WeightedEdge e);
        /**
         * Remove the given edge from the graph
         * @return true if the given edge has been removed, false otherwise (e.g. this edge does not exist)
         */
        virtual bool remove_edge(gfe::graph::Edge e);

        /**
         * Dump the content of the graph to given stream.
         */
        virtual void dump_ostream(std::ostream& out) const;

        /**
     * Perform a BFS from source_vertex_id to all the other vertices in the graph.
     * @param source_vertex_id the vertex where to start the search
     * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
     */
        virtual void bfs(std::uint64_t source_vertex_id, const char* dump2file = nullptr);

        /**
         * Execute the PageRank algorithm for the specified number of iterations.
         *
         * @param num_iterations the number of iterations to execute the algorithm
         * @param damping_factor weight for the PageRank algorithm, it affects the score associated to the sink nodes in the graphs
         * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
         */
        virtual void pagerank(std::uint64_t num_iterations, double damping_factor = 0.85, const char* dump2file = nullptr);

        /**
         * Weakly connected components (WCC), associate each node to a connected component of the graph
         * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
         */
        virtual void wcc(const char* dump2file = nullptr);

        /**
         * Community Detection using Label-Propagation. Associate a label to each vertex of the graph, according to its neighbours.
         * @param max_iterations max number of iterations to perform
         * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
         */
        virtual void cdlp(std::uint64_t max_iterations, const char* dump2file = nullptr);

        /**
         * Local clustering coefficient. Associate to each vertex the ratio between the number of its outgoing edges and the number of
         * possible remaining edges.
         * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
         */
        virtual void lcc(const char* dump2file = nullptr);

        /**
         * Single-source shortest paths. Compute the weight related to the shortest path from the source to any other vertex in the graph.
         * @param source_vertex_id the vertex where to start the search
         * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
         */
        virtual void sssp(std::uint64_t source_vertex_id, const char* dump2file = nullptr);

        virtual void generate_two_hops_neighbor_candidates(std::vector<std::uint64_t>&vertices);
        
        virtual void one_hop_neighbors(std::vector<std::uint64_t>&vertices);

        virtual void two_hop_neighbors(std::vector<std::uint64_t>&vertices);
            
    private:
        std::shared_ptr<pqxx::connection> get_conn() const;

    private:
        std::string _graph_name;
        const std::uint64_t INVALID_EXTID = 0;
        std::unique_ptr<ConnectionPool> _conn_pool;
    };

    class ConnectionPool {
        public:
        // 构造函数，初始化连接池
        ConnectionPool(const std::string& conninfo, std::size_t pool_size);

        ~ConnectionPool();

        // 获取连接
        std::shared_ptr<pqxx::connection> acquire();

        // 释放连接
        void release(std::shared_ptr<pqxx::connection> conn);

        // 调整连接池大小
        void resize(size_t size);

        private:
            std::unique_ptr<pqxx::connection> create_connection(const std::string& conninfo);

        private:
            std::string conninfo_;
            std::size_t pool_size_;
            std::queue<std::shared_ptr<pqxx::connection>> connections_;
            std::mutex mutex_;
            std::condition_variable condition_;
    };
}//namspace



