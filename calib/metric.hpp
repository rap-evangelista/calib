#include <vector>

namespace calib
{
    struct base_metric
    {
        public:
            virtual int value (int degree, unsigned long long int idx) = 0;
            virtual double _inner_prd_ (multivector& m1, multivector& m2) = 0;

            base_metric () {}
    };
}