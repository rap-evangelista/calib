#include <vector>

namespace calib
{
    struct base_metric
    {
        public:
            virtual void init () = 0;
            virtual std::vector <double> metric () = 0;
            virtual double _inner_prd_ (multivector& m1, multivector& m2) = 0;

            base_metric () {}
    };
}