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

            void add_elem (basis &e)
            {
                this-> elems. push_back (e);
            }
    };
}

#endif