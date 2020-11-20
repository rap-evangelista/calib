#ifndef CALIB_STRUCT_BASIS_OPERATORS
    #define CALIB_STRUCT_BASIS_OPERATORS

    #ifdef CALIB_MODE_HOST_DEVICE
        #include "basis_operators.cuh"
    #endif

namespace calib
{
    // # outer product possibilities for basis.

    const basis _outer_prd_ (basis& base1, basis& base2);   // # done.

    // # regressive product possibilities for basis.

    const basis _regr_prd_ (basis& base1, basis& base2);   // # done.

    // # inner product possibilities for basis.

    const basis _inner_prd_ (basis& base1, basis& base2, base_metric& metric);

    // # left contraction possibilities for basis.

    const basis _left_contr_ (basis& base1, basis& base2);

    // # individual operations possibilities for multivectors.

    const basis _rev_norm_ (basis& base);

    const basis _dual_ (basis& base);

    std::ostream& operator<< (std::ostream& os, basis& value)
    {
        os << value. direction () << " ";
        for (int i = 0; i < value. base_index. size (); i++)
        {
            os << "e" << value. base_index [i];
            if (i + 1 != value. base_index. size ()) os << "^";
        }

        return os;
    }

    // # host-only mode implementation.

    #ifdef CALIB_MODE_HOST_ONLY

        void cannonical_reordering (basis& base)
        {
            bool sorted = false;

            while (!sorted)
            {
                sorted = true;
                for (int i = 0; i < base. base_index. size () - 1; i++)
                {
                    // # se bases iguais, elas se eliminam
                    if (base. base_index [i] == base. base_index [i+1])
                    {
                        base. base_index. erase (base. base_index. begin () + i + 1);
                        base. base_index. erase (base. base_index. begin () + i --);
                        if (i < 0) break;
                        continue;
                    }

                    // # se base maior, troca as bases e inverte a orientação.
                    if (base. base_index [i] > base. base_index [i+1])
                    {
                        sorted = false;
                        int aux = base. base_index [i];
                        base. base_index [i] = base. base_index [i+1];
                        base. base_index [i+1] = aux;
                        base. orientation *= -1;
                    }
                }
            }
        }

    #endif

    // # host-device mode implementation.

    #ifdef CALIB_MODE_HOST_DEVICE
        void cannonical_reordering (basis& base)
        {
            // # verify if the size justify calling GPU or not.

            bridge_cannonical_reordering (base);
        }
    #endif

    // # always host-only implementation.

    template<typename T>
    const basis& operator* (T scalar, basis& base)
    {
        static_assert (std::is_arithmetic <T>::value, "An arithmetic type is required");
        base. magnitude = std::abs (scalar);
        base. orientation *= std::copysign (1, scalar);
        return base;
    }

    template<typename T>
    const basis& operator* (basis& base, T scalar)
    {
        static_assert (std::is_arithmetic <T>::value, "An arithmetic type is required");
        base. magnitude = std::abs (scalar);
        base. orientation *= std::copysign (1, scalar);
        return base;
    }

    const basis _outer_prd_ (basis& base1, basis& base2)
    {
        std::vector <int> base_index;
        base_index =  base1. base_index;

        for (int k = 0; k < base2. base_index. size (); k++)
        {
            base_index. emplace_back (base2. base_index [k]);
        }

        basis base3 = basis (base_index);
        base3. magnitude = base1. magnitude * base2. magnitude;
        base3. orientation = base1. orientation * base2. orientation;

        return base3;
    }

    const basis operator^ (basis& base1, basis& base2)
    {
        return _outer_prd_ (base1, base2);
    }

    const basis _regr_prd_ (basis& base1, basis& base2)
    {
        std::vector <int> base_index;

        int i = 0, j = 0;
        while (i < base1. base_index. size () && j < base2. base_index. size ())
        {
            if (base1. base_index [i] == base2. base_index [j])
            {
                base_index. push_back (base1. base_index [i]);
                i++; j++;
            } else {
                if (base1. base_index [i] < base2. base_index [j])
                {
                    i++;
                } else {
                    j++;
                }
            }
        }

        basis base3 = basis (base_index);
        base3. magnitude = base1. magnitude * base2. magnitude;
        base3. orientation = base1. orientation * base2. orientation;

        return base3;
    }

    const basis _inner_prd_ (basis& base1, basis& base2, base_metric& metric)
    {
        std::cout << "size of metric matrix: " << metric. matrix (). size () << std::endl;
        return basis ();
    }

    bool operator== (basis& base1, basis& base2)
    {
        return (base1. degree () == base2. degree ()) && (base1. unique_index () == base2. unique_index ());
    }
}

#endif