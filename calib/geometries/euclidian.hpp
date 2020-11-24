#ifndef CALIB_STRUCT_EUCLIDIAN_METRIC
    #define CALIB_STRUCT_EUCLIDIAN_METRIC

    #ifdef CALIB_MODE_HOST_DEVICE
        #include "euclidian.cuh"
    #endif

namespace calib
{
    struct euclidian_metric : public base_metric
    {
        public:
            std::vector <double> metric_mtx;

            virtual void init () override
            {
                this-> metric_mtx = std::vector <double> ({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
            }

            virtual std::vector <double> metric () override
            {
                return this-> metric_mtx;
            }

        #ifdef CALIB_MODE_HOST_ONLY
            virtual double _inner_prd_ (multivector& m1, multivector& m2) override
            {
                return 0;
            }
        #endif

        #ifdef CALIB_MODE_HOST_DEVICE
            virtual double _inner_prd_ (multivector& m1, multivector& m2) override
            {
                return bridge_metric_inner_prd (m1, m2, this-> metric_mtx);
            }
        #endif



            euclidian_metric () : base_metric () {init ();}
    };
}

#endif