#include <cmath>
#include <type_traits>

namespace calib
{
    template<typename T>
    basis& operator* (T scalar, basis& base)
    {
        static_assert (std::is_arithmetic <T>::value, "An arithmetic type is required");
        base. magnitude = std::abs (scalar);
        base. orientation = std::copysign (1, scalar);
        return base;
    }

    template<typename T>
    basis& operator* (basis& base, T scalar)
    {
        static_assert (std::is_arithmetic <T>::value, "An arithmetic type is required");
        base. magnitude = std::abs (scalar);
        base. orientation = std::copysign (1, scalar);
        return base;
    }

    const basis operator^ (basis& base1, basis& base2)
    {
        basis base3 = basis ();
        base3. magnitude = base1. magnitude * base2. magnitude;
        base3. orientation = base1. orientation * base2. orientation;

        // # une os index e rearruma.
        base3. base_index = base1. base_index;

        return base3;
    }

    multivector& operator^ (multivector& v1, multivector& v2)
    {
        return multivector ();
    }

    std::ostream& operator<<(std::ostream& os, basis& value)
    {
        os << value. direction () << " ";
        for (int i = 0; i < value. base_index. size (); i++)
        {
            os << "e" << value. base_index [i];
            if (i + 1 != value. base_index. size ()) os << "^";
        }

        return os;
    }
}