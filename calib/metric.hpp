#include <vector>

namespace calib
{
    void metric_init (double * metric_matrix);

    struct metric
    {
        public:
            //static bool initialized;
            //static double metric_matrix [];

            metric ()
            {
                //if (!initialized) metric_init (metric_matrix);
                //initialized = true;
            }
    };
}