namespace calib
{
    struct euclidian_metric : public base_metric
    {
        public:
            virtual void init () override
            {
                //metric_matrix [0] = 1;
                std::cout << "euclidian_metric" << std::endl;
            }

            virtual std::vector <double> matrix () override
            {
                std::vector <double> mtx ({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
                return mtx;
            }

            euclidian_metric () : base_metric () {init ();}
    };
}