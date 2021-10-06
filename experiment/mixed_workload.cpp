//
// Created by per on 29.09.21.
//

#include "mixed_workload.hpp"

#include <future>
#include <chrono>

#include "graphalytics.hpp"
#include "aging2_experiment.hpp"
#include "mixed_workload_result.hpp"

namespace gfe::experiment {

    using namespace std;

    MixedWorkloadResult MixedWorkload::execute() {
      auto aging_result_future = std::async(std::launch::async, &Aging2Experiment::execute, &m_aging_experiment);

      chrono::seconds progress_check_interval( 1 );
      this_thread::sleep_for( progress_check_interval );  // Poor mans synchronization to ensure AgingExperiment was able to setup the master etc
      while (true) {
        if (m_aging_experiment.progress_so_far() > 0.1) { // The graph reached its final size
          break;
        }
        this_thread::sleep_for( progress_check_interval ) ;
      }
      cout << "Graph reached final size, executing graphaltyics now" << endl;

      // TODO change this to also work for LCC, generally this solution is very ugly
#if defined(HAVE_OPENMP)
      if(configuration().num_threads(ThreadsType::THREADS_READ) != 0 ){
                        LOG("[driver] OpenMP, number of threads for the Graphalytics suite: " << configuration().num_threads(ThreadsType::THREADS_READ));
                        omp_set_num_threads(configuration().num_threads(ThreadsType::THREADS_READ));
                    }
#endif

      m_graphalytics.execute();

      if (m_aging_experiment.progress_so_far() > 0.98) {
        cerr << "Analytics stopped after updates potentially stopped. This run is invalid." << endl;
      }
      aging_result_future.wait();
      auto aging_result = aging_result_future.get();

      return MixedWorkloadResult { aging_result, m_graphalytics };
    }

}