#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <iostream>
#include <random>
#include <vector>
#include <cstdarg>

#define DEFAULT_SPACE_DIM 30

// # define running mode before the 'core'.

#include "calib/mode/host_device.hpp"
//#include "calib/mode/host_only.hpp"

// # choose default metric space.

//#include "calib/geometries/euclidian.hpp"

#include "calib/core.hpp"

using namespace calib;

multivector random_generate_multivector ()
{
    std::srand (time (0));

    multivector m1;

    // # random number of elements .

    int dim = DEFAULT_SPACE_DIM;

    for (int i = 0; i < dim; i++)
    {
        // # random degree of current element.

        int degree = (std::rand() % (dim))+1;

        std::vector <int> indices;

        for (int j = 0; j < degree; j++)
        {
            int random = (std::rand() % (DEFAULT_SPACE_DIM))+1;
            indices. push_back (random != 0 ? random : 1);
        }

        basis base (indices);
        base. magnitude = std::rand () % 100;
        base. orientation = std::copysign (1, (std::rand() % 2) - 1);

        m1. add_elem (base);
    }

    cannonical_reordering (m1);

    return m1;
}

int main (int argv, char * argc)
{
    auto e1 = basis ({2,4,3}) * -3;
    auto e2 = basis ({2}) * 2;
    auto e3 = basis ({3,2,1}) * -7;

    //std::cout << "e1: " << e1 << std::endl;
    //std::cout << "e2: " << e2 << std::endl;

    multivector m1 = e1 + e2;
    multivector m2 = e2 + e3;

    auto m3 = _outer_prd_ (m1, m2);
    auto m4 = _regr_prd_ (m1, m2);
    float inner_prd = _inner_prd_ (m1, m2, euclidian_metric ());

    //std::cout << "m1: " << m1 << std::endl;
    //std::cout << "m2: " << m2 << std::endl;

    //std::cout << "m1 . m2: " << inner_prd << std::endl;

#ifdef CALIB_MODE_HOST_DEVICE
    // # init constant memory

    float * fact;
    fact = (float *) std::malloc (CALIB_FACTORIAL_TABLE_SIZE * sizeof (float));

    fact [0] = 1;
    for (int i = 1; i < CALIB_FACTORIAL_TABLE_SIZE; i++)
    {
        fact [i] = i * fact [i-1];
    }

    cudaMemcpyToSymbol (fact_, &fact, CALIB_FACTORIAL_TABLE_SIZE * sizeof(float));
#endif

    // # rules of the experiment.

    int n = 1000;
    double mean_time = 0;


    for (int i = 0; i < n; i++)
    {
        // # random generation of

        multivector m1 = random_generate_multivector ();
        multivector m2 = random_generate_multivector ();

        #ifdef CALIB_MODE_HOST_ONLY
            std::clock_t c_start = std::clock();
        #endif

        #ifdef CALIB_MODE_HOST_DEVICE
            cudaEvent_t start, stop;

            cudaEventCreate (&start);
            cudaEventCreate (&stop);

            cudaEventRecord (start, 0);
        #endif

        //m1 ^ m2;
        //_regr_prd_ (m1, m2);
        _inner_prd_ (m1, m2, euclidian_metric ());

        #ifdef CALIB_MODE_HOST_ONLY
            std::clock_t c_end   = std::clock();
            mean_time += 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
        #endif

        #ifdef CALIB_MODE_HOST_DEVICE
            cudaEventRecord (stop, 0);

            cudaEventSynchronize (stop);

            float time_elapsed_ms;
            cudaEventElapsedTime (&time_elapsed_ms, start, stop);
            mean_time += time_elapsed_ms;
            
            cudaEventDestroy (start);
            cudaEventDestroy (stop);
        #endif

        std::cout << ".";
    }

#ifdef CALIB_MODE_HOST_DEVICE
    cudaFree (&fact_);
#endif

    std::cout << std::endl;

    std::cout << "Experiment with " << DEFAULT_SPACE_DIM << " dimensions." << std::endl;
    std::cout << "> running after " << n << " times." << std::endl;
    std::cout << "> mean time: " << mean_time / n << " ms." << std::endl;
}