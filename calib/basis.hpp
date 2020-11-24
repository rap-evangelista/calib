#include <cstdarg>
#include <vector>

#ifndef CALIB_STRUCT_BASIS
    #define CALIB_STRUCT_BASIS

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

            int unique_index ()
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

        private:
            int cantor_pairing (int x, int y)
            {
                return (x * x + 3 * x + 2 * x * y + y + y * y) / 2;
            }
    };
}

#endif