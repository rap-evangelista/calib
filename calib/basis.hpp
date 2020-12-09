#include <cstdarg>
#include <vector>

#ifndef CALIB_STRUCT_BASIS
    #define CALIB_STRUCT_BASIS

#include <algorithm>

namespace calib
{
    struct basis;
    void cannonical_reordering (basis& base);

    struct basis
    {
        public:
            std::vector <int> base_index;   // # a list of indices representing the base vectors (e1 -> {1}, e2 -> {2}, e1^e2 -> {1,2}, ...).
            int orientation = 1;            // # the way direction how a particle moves through a support.
            double magnitude = 1;           // # the speed of the particle.

            basis () {}
            basis (std::vector <int> args) {this-> init (args);}

            void init (std::vector <int> args)
            {
                this-> base_index = args;

                // # necessário verificar se algum dos argumentos é maior que a dimensão default do espaço.
                for (auto v : this-> base_index)
                    if (v > DEFAULT_SPACE_DIM)
                        std::cout << "A basis can't have an index greater than the default space dimension (which is " << DEFAULT_SPACE_DIM << ").";

                // # necessário verificar ordem dos argumentos.
                // # se não estiver em ordem crescente, aplica a reordenação canônica.
                // # se houverem bases repetidas, o blade cairá de grau.
                cannonical_reordering (*this);
            }

            int degree ()
            {
                return this-> base_index. size ();
            }

            double direction ()
            {
                return this-> orientation * this-> magnitude;
            }
            
            unsigned long long int cantor_pairing (unsigned long long int x, unsigned long long int y)
            {
                return (x * x + 3 * x + 2 * x * y + y + y * y) / 2;
            }

            void recover_from_degree_and_unique_index (int degree, unsigned long long int unique_index)
            {
                this-> base_index. clear ();
                
                if (degree == 1)
                    this-> base_index. push_back ((int) unique_index);
                
                if (degree > 1)
                {
                    unsigned long long int next_unique_index = unique_index;

                    for (int i = 0; i < degree; i++)
                    {
                        if (i == degree - 1)
                        {
                            this-> base_index. push_back ((int) next_unique_index);
                            break;
                        }

                        unsigned long long int aux = std::floor ((-1 + std::sqrt (1 + 8 * next_unique_index)) / 2.0f);

                        unsigned long long int x = next_unique_index - (aux * (1 + aux) / 2.0f);
                        unsigned long long int y = (aux * (3 + aux) / 2) - next_unique_index;

                        next_unique_index = x;

                        this-> base_index. push_back ((int) y);
                    }

                    std::reverse (std::begin (this-> base_index), std::end (this-> base_index));
                }
            }

            unsigned long long int unique_index ()
            {
                if (this-> base_index. size () == 0)
                    return 0;

                int value = this-> base_index [0];

                if (this-> base_index. size () == 1)
                    return value;
                
                for (int i = 1; i < this-> base_index. size (); i++)
                {
                    value = this-> cantor_pairing (value, this-> base_index [i]);
                }

                return value;
            }
    };
}

#endif