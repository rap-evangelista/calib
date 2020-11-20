#include <vector>

namespace calib
{
    struct base_metric
    {
        public:
            virtual void init () = 0;
            virtual std::vector <double> matrix () = 0;

            base_metric () {}
    };
}