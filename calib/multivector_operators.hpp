#ifndef CALIB_STRUCT_MULTIVECTOR_OPERATORS
    #define CALIB_STRUCT_MULTIVECTOR_OPERATORS

    #ifdef CALIB_MODE_HOST_DEVICE
        #include "multivector_operators.cuh"
    #endif

namespace calib
{
    // # scalar multiplication possibilities for multivectors.

    template<typename T>
    multivector& operator* (T scalar, multivector& mv);

    template<typename T>
    multivector& operator* (multivector& mv, T scalar);

    // # summation possibilities between multivectors and basis.

    const multivector operator+ (basis& base1, basis& base2);

    const multivector operator+ (multivector& mv, basis& base);

    const multivector operator+ (basis& base, multivector& mv);

    const multivector operator+ (multivector& mv1, multivector& mv2);

    template<typename T>
    const basis operator+ (T scalar, basis& base);

    template<typename T>
    const basis operator+ (basis& base, T scalar);

    template<typename T>
    const multivector operator+ (T scalar, multivector& mv);

    template<typename T>
    const multivector operator+ (multivector& mv, T scalar);

    // # outer product possibilities between multivectors and basis.

    const multivector& operator^ (basis& base, multivector& mv);

    const multivector& operator^ (multivector& mv, basis& base);

    const multivector& operator^ (multivector& v1, multivector& v2);

    const multivector _outer_prd_ (basis& base, multivector& mv);

    const multivector _outer_prd_ (multivector& mv, basis& base);

    const multivector _outer_prd_ (multivector& v1, multivector& v2);

    // # regressive product possibilities between multivectors and basis.

    const multivector _regr_prd_ (basis& base, multivector& mv);

    const multivector _regr_prd_ (multivector& mv, basis& base);

    const multivector _regr_prd_ (multivector& v1, multivector& v2);

    // # inner product possibilities between multivectors and basis.

    const multivector _inner_prd_ (basis& base, multivector& mv);

    const multivector _inner_prd_ (multivector& mv, basis& base);

    const multivector _inner_prd_ (multivector& v1, multivector& v2);

    // # left contraction possibilities between multivectors and basis.

    const multivector _left_contr_ (basis& base, multivector& mv);

    const multivector _left_contr_ (multivector& mv, basis& base);

    const multivector _left_contr_ (multivector& v1, multivector& v2);

    // # geometric product possibilities between multivectors and basis.

    const multivector _geom_prd_ (basis& base, multivector& mv);

    const multivector _geom_prd_ (multivector& mv, basis& base);

    const multivector _geom_prd_ (multivector& v1, multivector& v2);

    // # individual operations possibilities for multivectors.

    const multivector _rev_norm_ (multivector& mv);

    const multivector _dual_ (multivector& mv);

    std::ostream& operator<< (std::ostream& os, multivector& value)
    {
        for (int i = 0; i < value. elems. size (); i++)
        {
            os << "(" << value. elems [i] << ")";
            if (i + 1 != value. elems. size ()) os << "+";
        }

        return os;
    }

    // # host-only mode implementation.

    #ifdef CALIB_MODE_HOST_ONLY

        const multivector operator+ (basis& base1, basis& base2)
        {
            multivector mv = multivector ();

            mv. add_elem (base1);
            mv. add_elem (base2);

            // # rearrange.

            for (int i = 0; i < mv. elems. size () - 1; i++)
            {
                for (int k = (i+1); k < mv. elems. size (); k++)
                {
                    if (mv. elems [k] == mv. elems [i])
                    {
                        float scalar = mv. elems [i]. magnitude + mv. elems [k]. direction ();
                        mv. elems [i]. magnitude = std::abs (scalar);
                        mv. elems [i]. orientation *= std::copysign (1, scalar);
                        mv. elems. erase (mv. elems. begin () + k --);
                    }
                }
            }

            return mv;
        }

        // # outer product.

        const multivector _outer_prd_ (multivector& v1, multivector& v2)
        {
            multivector v3 = multivector ();

            for (auto base1 : v1. elems)
            {
                for (auto base2 : v2. elems)
                {
                    basis base_ = base1 ^ base2;
                    v3. add_elem (base_);
                }
            }

            // # reduct.

            return v3;
        }

        const multivector& operator^ (multivector& v1, multivector& v2)
        {
            return _outer_prd_ (v1, v2);
        }

        // # regressive product.

        const multivector _regr_prd_ (multivector& v1, multivector& v2)
        {
            multivector v3 = multivector ();

            for (auto base1 : v1. elems)
            {
                for (auto base2 : v2. elems)
                {
                    basis base_ = _regr_prd_ (base1, base2);
                    v3. add_elem (base_);
                }
            }

            // # reduct.

            return v3;
        }

        // # inner product.

    #endif

    // # host-device mode implementation.

    #ifdef CALIB_MODE_HOST_DEVICE

        const multivector operator+ (basis& base1, basis& base2)
        {
            return multivector ();
        }

        const multivector operator+ (multivector& mv1, multivector& mv2)
        {
            return bridge_sum_multivectors (mv1, mv2);
            //return multivector ();
        }

        // # external product.

        const multivector& operator^ (multivector& mv1, multivector& mv2)
        {
            return bridge_outer_prd_multivectors (mv1, mv2);
            //return multivector ();
        }

        // # inner product.

    #endif

    // # always host-only implementation.

}

#endif