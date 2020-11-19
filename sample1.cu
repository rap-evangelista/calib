#include <stdio.h>
#include <iostream>

#define DEFAULT_SPACE_DIM 3

// # define running mode before the 'core'.

//#include "calib/mode/host_device.hpp"
#include "calib/mode/host_only.hpp"

#include "calib/core.hpp"

using namespace calib;

int main (int argv, char * argc)
{
    auto e1 = basis ({3,2,1}) * -3;
    auto e2 = -2 * basis ({2,1,3,1});
    auto e3 = e1 ^ e2;

    std::cout << e1 << std::endl;
    std::cout << e2 << std::endl;
    std::cout << e3 << std::endl;

    auto m1 = e1 + e3;

    std::cout << m1 << std::endl;
}