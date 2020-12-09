#ifndef CALIB_STRUCT_MULTIVECTOR
    #define CALIB_STRUCT_MULTIVECTOR

namespace calib
{
    struct multivector;
    void cannonical_reordering (multivector& base);
    
    struct multivector
    {
        public:
            std::vector <basis> elems;

            multivector () {}

            void clear ()
            {
                elems. clear ();
            }

            void add_elem (basis &e)
            {
                this-> elems. push_back (e);
            }

            basis search (int degree, unsigned long long int index)
            {
                for (auto base : this-> elems)
                {
                    if (base. degree () == degree && base. unique_index () == index) return base;
                }

                basis null;
                null. magnitude = 0;
                null. orientation = -1;
                return null;
            }
    };
}

#endif