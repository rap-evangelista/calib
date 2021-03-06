#ifndef CALIB_STRUCT_BASIS_OPERATORS
    #define CALIB_STRUCT_BASIS_OPERATORS

    #ifdef CALIB_MODE_HOST_DEVICE
        #include "basis_operators.cuh"
    #endif

namespace calib
{
    // # outer product possibilities for basis.

    basis _outer_prd_ (basis& base1, basis& base2);   // # done.

    // # regressive product possibilities for basis.

    basis _regr_prd_ (basis& base1, basis& base2);   // # done.

    // # individual operations possibilities for multivectors.

    basis _reverse_ (basis const base);    // # done.

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

            if (base. base_index. size () <= 0) return;

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
    basis operator* (T scalar, basis& base)
    {
        static_assert (std::is_arithmetic <T>::value, "An arithmetic type is required");
        base. magnitude = std::abs (scalar);
        base. orientation *= std::copysign (1, scalar);
        return base;
    }

    template<typename T>
    basis operator* (basis& base, T scalar)
    {
        static_assert (std::is_arithmetic <T>::value, "An arithmetic type is required");
        base. magnitude = std::abs (scalar);
        base. orientation *= std::copysign (1, scalar);
        return base;
    }

    basis _outer_prd_ (basis& base1, basis& base2)
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

    basis operator^ (basis& base1, basis& base2)
    {
        return _outer_prd_ (base1, base2);
    }

    basis _regr_prd_ (basis& base1, basis& base2)
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

    basis _reverse_ (basis const base_)
    {
        std::vector <int> base_index;

        for (int i = base_. base_index. size () - 1; i >= 0; i++)
        {
            base_index. emplace_back (base_. base_index [i]);
        }

        basis base = basis (base_index);
        base. magnitude = base_. magnitude;
        base. orientation = pow (-1, DEFAULT_SPACE_DIM * (DEFAULT_SPACE_DIM - 1) / 2) * base_. orientation;

        return base;
    }

    bool operator== (basis& base1, basis& base2)
    {
        return (base1. degree () == base2. degree ()) && (base1. unique_index () == base2. unique_index ());
    }

    bool operator> (basis& base1, basis& base2)
    {
        return (base1. degree () > base2. degree ()) || (base1. degree () == base2. degree ()) && (base1. unique_index () > base2. unique_index ());
    }

    bool operator< (basis& base1, basis& base2)
    {
        return (base1. degree () < base2. degree ()) || (base1. degree () == base2. degree ()) && (base1. unique_index () < base2. unique_index ());
    }
}

#endif