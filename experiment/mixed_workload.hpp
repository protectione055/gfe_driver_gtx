#ifndef GFE_DRIVER_MIXED_WORKLOAD_H
#define GFE_DRIVER_MIXED_WORKLOAD_H

namespace gfe::experiment { class Aging2Experiment; }
namespace gfe::experiment { class GraphalyticsSequential; }
namespace gfe::experiment { class MixedWorkloadResult; }

namespace gfe::experiment {

    class MixedWorkload {
    public:
        MixedWorkload(Aging2Experiment& aging_experiment, GraphalyticsSequential& graphalytics)
          : m_aging_experiment(aging_experiment), m_graphalytics(graphalytics) {}

        MixedWorkloadResult execute();
    private:
        Aging2Experiment& m_aging_experiment;
        GraphalyticsSequential& m_graphalytics;
    };

}

#endif //GFE_DRIVER_MIXED_WORKLOAD_H
