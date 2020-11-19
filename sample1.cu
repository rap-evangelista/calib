#include <stdio.h>
#include <iostream>

#define DEFAULT_SPACE_DIM 3

// # define running mode before the 'core'.

#include "calib/mode/host_device.hpp"
//#include "calib/mode/host_only.hpp"
//    or "calib/mode/host_only.hpp"

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

    //auto m1 = e1 + e3;
    //auto m2 = e1 + e2;
    //auto m3 = m1 ^ m2;
    //mv. elements [0] = e (1);
    //auto x = e (1);
    //std::cout << m1 << std::endl;
    //std::cout << m2 << std::endl;
    //std::cout << m3 << std::endl;
}