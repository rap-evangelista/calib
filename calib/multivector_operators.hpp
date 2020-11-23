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

    const multivector operator+ (basis& base1, basis& base2);   // # done.

    const multivector operator+ (multivector& mv, basis& base); // # done.

    const multivector operator+ (basis& base, multivector& mv); // # done.

    const multivector operator+ (multivector& mv1, multivector& mv2);   // # done.

    template<typename T>
    const basis operator+ (T scalar, basis& base);  // # done.

    template<typename T>
    const basis operator+ (basis& base, T scalar);  // # done.

    template<typename T>
    const multivector operator+ (T scalar, multivector& mv);    // # done.

    template<typename T>
    const multivector operator+ (multivector& mv, T scalar);    // # done.

    // # outer product possibilities between multivectors and basis.

    const multivector& operator^ (basis& base, multivector& mv);    // # done.

    const multivector& operator^ (multivector& mv, basis& base);    // # done.

    const multivector& operator^ (multivector& v1, multivector& v2);   // # done.

    const multivector _outer_prd_ (basis& base, multivector& mv);    // # done.

    const multivector _outer_prd_ (multivector& mv, basis& base);    // # done.

    const multivector _outer_prd_ (multivector& v1, multivector& v2);   // # done.

    // # regressive product possibilities between multivectors and basis.

    const multivector _regr_prd_ (basis& base, multivector& mv);    // # done.

    const multivector _regr_prd_ (multivector& mv, basis& base);    // # done.

    const multivector _regr_prd_ (multivector& v1, multivector& v2);    // # done.

    // # inner product possibilities between multivectors and basis.

    const double _inner_prd_ (basis& base1, basis& base2, base_metric& metric);

    const double _inner_prd_ (basis& base, multivector& mv, base_metric& metric);

    const double _inner_prd_ (multivector& mv, basis& base, base_metric& metric);

    const double _inner_prd_ (multivector& v1, multivector& v2, base_metric& metric);

    // # left contraction possibilities between multivectors and basis.

    const multivector _left_contr_ (basis& base1, basis& base2);

    const multivector _left_contr_ (basis& base, multivector& mv, base_metric& metric);

    const multivector _left_contr_ (multivector& mv, basis& base, base_metric& metric);

    const multivector _left_contr_ (multivector& v1, multivector& v2, base_metric& metric);

    // # geometric product possibilities between multivectors and basis.

    const multivector _geom_prd_ (basis& base1, basis& base2, base_metric& metric);

    const multivector _geom_prd_ (basis& base, multivector& mv, base_metric& metric);

    const multivector _geom_prd_ (multivector& mv, basis& base, base_metric& metric);

    const multivector _geom_prd_ (multivector& v1, multivector& v2, base_metric& metric);

    // # individual operations possibilities for multivectors.

    const multivector _rev_norm_ (multivector& mv);

    const multivector _dual_ (multivector& mv);

    // # host-only mode implementation.

    #ifdef CALIB_MODE_HOST_ONLY

        void cannonical_reordering (multivector& mv)
        {
            bool sorted = false;

            while (!sorted)
            {
                sorted = true;
                for (int i = 0; i < mv. elems. size () - 1; i++)
                {
                    // # se bases iguais, elas se eliminam
                    if (mv. elems [i] == mv. elems [i+1])
                    {
                        // # merge equal basis.
                        auto scalar = mv. elems [i]. direction () + mv. elems [i+1]. direction ();
                        mv. elems [i]. magnitude = std::abs (scalar);
                        mv. elems [i]. orientation *= std::copysign (1, scalar);

                        mv. elems. erase (mv. elems. begin () + i + 1);

                        continue;
                    }

                    // # se base maior, troca as bases e inverte a orientação.
                    if (mv. elems [i] > mv. elems [i+1])
                    {
                        sorted = false;
                        basis aux = mv. elems [i];
                        mv. elems [i] = mv. elems [i+1];
                        mv. elems [i+1] = aux;
                    }
                }
            }
        }

        const multivector operator+ (basis& base1, basis& base2)
        {
            multivector mv = multivector ();

            mv. add_elem (base1);
            mv. add_elem (base2);

            // # rearrange.
            cannonical_reordering (mv);

            return mv;
        }

        const multivector operator+ (multivector& mv_, basis& base)
        {
            multivector mv = multivector ();

            for (auto base : mv_. elems)
            {
                mv. elems. emplace_back (base);
            }

            mv. add_elem (base);

            cannonical_reordering (mv);

            return mv;
        }

        const multivector operator+ (basis& base, multivector& mv_)
        {
            return (mv_ + base);
        }

        const multivector operator+ (multivector& mv1, multivector& mv2)
        {
            multivector mv = multivector ();

            for (auto base : mv1. elems)
            {
                mv. elems. emplace_back (base);
            }

            for (auto base : mv2. elems)
            {
                mv. elems. emplace_back (base);
            }

            cannonical_reordering (mv);

            return mv;
        }

        template<typename T>
        const multivector operator+ (T scalar, basis& base)
        {
            basis basex = basis ();

            basex. magnitude = std::abs (scalar);
            basex. orientation *= std::copysign (1, scalar);

            return (basex + base);
        }

        template<typename T>
        const multivector operator+ (basis& base, T scalar)
        {
            return (scalar + base);
        }

        template<typename T>
        const multivector operator+ (T scalar, multivector& mv)
        {
            basis basex = basis ();

            basex. magnitude = std::abs (scalar);
            basex. orientation *= std::copysign (1, scalar);

            return (basex + mv);
        }

        template<typename T>
        const multivector operator+ (multivector& mv, T scalar)
        {
            return (scalar + mv);
        }

        // # outer product.

        const multivector _outer_prd_ (basis& base, multivector& mv_)
        {
            multivector mv = multivector ();

            mv. add_elem (base);

            for (basis base_ : mv_. elems)
            {
                mv. add_elem (base_);
            }

            cannonical_reordering (mv);

            return mv;
        }

        const multivector _outer_prd_ (multivector& mv, basis& base)
        {
            return _outer_prd_ (base, mv);
        }

        const multivector _outer_prd_ (multivector& mv1, multivector& mv2)
        {
            multivector mv = multivector ();

            for (auto base1 : mv1. elems)
            {
                for (auto base2 : mv2. elems)
                {
                    basis base_ = base1 ^ base2;
                    mv. add_elem (base_);
                }
            }

            cannonical_reordering (mv);

            return mv;
        }

        const multivector& operator^ (basis& base, multivector& mv)
        {
            return _outer_prd_ (base, mv);
        }

        const multivector& operator^ (multivector& mv, basis& base)
        {
            return _outer_prd_ (base, mv);
        }

        const multivector& operator^ (multivector& mv1, multivector& mv2)
        {
            return _outer_prd_ (mv1, mv2);
        }

        // # regressive product.

        const multivector _regr_prd_ (basis& base, multivector& mv_)
        {
            multivector mv = multivector ();

            for (auto base1 : mv_. elems)
            {
                basis base_ = _regr_prd_ (base1, base);
                mv. add_elem (base_);
            }

            cannonical_reordering (mv);

            return mv;
        }

        const multivector _regr_prd_ (multivector& mv, basis& base)
        {
            return _regr_prd_ (base, mv);
        }

        const multivector _regr_prd_ (multivector& mv1, multivector& mv2)
        {
            multivector mv = multivector ();

            for (auto base1 : mv1. elems)
            {
                for (auto base2 : mv2. elems)
                {
                    basis base_ = _regr_prd_ (base1, base2);
                    mv. add_elem (base_);
                }
            }

            cannonical_reordering (mv);

            return mv;
        }

        const double _inner_prd_ (basis& base1, basis& base2, base_metric& metric)
        {
            multivector m1 = multivector ();
            multivector m2 = multivector ();

            m1. add_elem (base1);
            m2. add_elem (base2);

            return metric. _inner_prd_ (m1, m2);
        }

        // # left contraction.

        const multivector _left_contr_ (basis& base1, basis& base2)
        {
            //basis dual2 = _dual_ (base2);
            //basis base1_dual2 = (base1 ^ dual2);
            //basis udual = _undual_ (base1_dual2);

            multivector m1 = multivector ();
            //m1. add_elem (udual);

            return m1;
        }

        // # geometric product.

        const multivector _geom_prd_ (basis& base1, basis& base2, base_metric& metric)
        {
            double scalar = _inner_prd_ (base1, base2, metric);
            basis base_r = _outer_prd_ (base1, base2);

            multivector m1 = multivector ();
            m1. add_elem (base_r);

            return m1;
        }

    #endif

    // # host-device mode implementation.

    #ifdef CALIB_MODE_HOST_DEVICE

        void cannonical_reordering (multivector& mv)
        {
            
        }

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

    std::ostream& operator<< (std::ostream& os, multivector& value)
    {
        for (int i = 0; i < value. elems. size (); i++)
        {
            os << "(" << value. elems [i] << ")";
            if (i + 1 != value. elems. size ()) os << "+";
        }

        return os;
    }

}

#endif