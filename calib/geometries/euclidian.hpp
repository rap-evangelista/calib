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
            int p = 5;      // # 0 <= p <= q
            int q = 8;      // # 0 <= q <= 2^n

            float fact (int x)
            {
                return (x <= 1) ? 1 : x * fact (x-1);
            }

            void find_in_combinations (int k, int ini, int fim, int p, int q, unsigned long long int idx, int start_index, int & relative_pos, int & global_pos, std::vector <int> & all_indices, std::vector <int> & aux_combinations, int & result)
            {
                if (k == 0)
                {
                    unsigned long long int index = aux_combinations [0];

                    for (int i = 1; i < aux_combinations. size (); i++)
                    {
                        index = basis (). cantor_pairing (index, aux_combinations [i]);
                    }

                    if (index == idx)
                    {
                        global_pos = (int) ((float) relative_pos);
                        if (global_pos <= p) result = 1;
                        if (global_pos > q) result = 0;
                        if (global_pos > p && global_pos <= q) result = -1;
                    }

                    relative_pos += 1;

                    return;
                }

                for (int i = start_index; i <= all_indices. size () - k; i++)
                {
                    aux_combinations. push_back (all_indices [i]);
                    find_in_combinations (k - 1, ini, fim, p, q, idx, i + 1, relative_pos, global_pos, all_indices, aux_combinations, result);
                    aux_combinations. pop_back ();
                }
            }

            int find_in_combinations (int degree, int ini, int fim, int p, int q, unsigned long long int idx)
            {
                int relative_pos = ini;
                int global_pos = ini + 1;
                int result = 0;

                std::vector <int> all_indices, aux_combinations;
                for (int i = 1; i <= DEFAULT_SPACE_DIM; i++)
                {
                    all_indices. push_back (i);
                }

                find_in_combinations (degree, ini, fim, p, q, idx, 0, relative_pos, global_pos, all_indices, aux_combinations, result);

                return result;
            }

            virtual int value (int degree, unsigned long long int idx) override
            {
                int ini = 0;
                int fim = 0;

                for (int k = 0; k <= degree - 1; k++)
                {
                    ini += fact (DEFAULT_SPACE_DIM) / (fact (DEFAULT_SPACE_DIM - k) * fact (k));
                }

                if (ini >= this-> q) return 0;   // # o intervalo inicia-se após q.

                fim = ini + fact (DEFAULT_SPACE_DIM) / (fact (DEFAULT_SPACE_DIM - degree) * fact (degree));

                if (fim < this-> p) return 1;   // # o intervalo finaliza-se antes de p.

                // # o intervalo cruza (p, q).

                return find_in_combinations (degree, ini, fim, p, q, idx);
            }

        #ifdef CALIB_MODE_HOST_ONLY
            virtual double _inner_prd_ (multivector& m1, multivector& m2) override
            {
                double result = 0;

                for (auto base1 : m1. elems)
                {
                    auto base2 = m2. search (base1. degree (), base1. unique_index ());

                    if (base2. direction () != 0)
                    {
                        result += base1. direction () * value (base1. degree (), base1. unique_index ()) * base2. direction ();
                    }
                }

                return result;
            }
        #endif

        #ifdef CALIB_MODE_HOST_DEVICE
            virtual double _inner_prd_ (multivector& m1, multivector& m2) override
            {
                return bridge_metric_inner_prd (m1, m2);
            }
        #endif

            euclidian_metric () : base_metric () {}
    };
}

#endif