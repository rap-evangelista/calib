namespace calib
{
    struct homogeneous_metric : public base_metric
    {
        public:
            virtual void init () override
            {
                //metric_matrix [0] = 1;
                std::cout << "homogeneous_metric" << std::endl;
            }

            virtual std::vector <double> matrix () override
            {
                std::vector <double> mtx ({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
                return mtx;
            }

            virtual double _inner_prd_ (multivector& m1, multivector& m2)
            {
                return 0;
            }

            homogeneous_metric () : base_metric () {init ();}
    };
}