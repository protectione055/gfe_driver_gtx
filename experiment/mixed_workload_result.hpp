#ifndef GFE_DRIVER_MIXED_WORKLOAD_RESULT_H
#define GFE_DRIVER_MIXED_WORKLOAD_RESULT_H

namespace gfe::experiment { class GraphalyticsSequential; }
namespace gfe::experiment { class Aging2Result; }
namespace common { class Database; }

namespace gfe::experiment {

    class MixedWorkloadResult {
    public:
        MixedWorkloadResult(Aging2Result& aging_result, GraphalyticsSequential& analytics);

        void save(common::Database* db);

    private:
        Aging2Result& m_aging_result;
        GraphalyticsSequential& m_graphalytics;
    };

}

#endif //GFE_DRIVER_MIXED_WORKLOAD_RESULT_H
