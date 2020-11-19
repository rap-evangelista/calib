# CALib - A Simple Clifford Algebra Library
CALib is a C++ library to manipulate subspace as primitives, based on Clifford Algebra. The code structure is hardly following the 'keep it simple' policy.

## Prime Goal
This library is still under construction and will be presented as a final project for the 'GPU Architecture and Programming' discipline, under the watch of the professor D.Sc. Esteban Clua, at UFF.

## Available Functionalities
 - Manipulate *subspaces* in a fixed *n*-dimensional space.
 - Automatic **cannonical reordering** * to mantain an organized basis.
 - **Outer product** *.
 - **Regressive product** *.
 - **Inner product** * with a customizable metric matrix (Euclidian by default).
 - **Left contraction** * operator.
 - **Geometric Product** *.
 - **Reverse norm** operator.
 - **Dual** operator.
 
 '* CUDA optimized.

## Code Organization
The code is organized as follows:

 1. Any file `sampleX.cu` contains the main function and provides an example of how to use the library.
 2. The subdirectory `calib/*` contains all the files used by the library.
 3. The subdirectory `calib/mode/*` contains files that indicate the possible ways of execution ('til now, Host-Only, for CPU version, and Host-Device, for hybrid CPU-GPU version).
 4. There are two structs to be used `basis.hpp` and `multivector.hpp`.
 5. Each one of them has a separated file containing its respectives `x_operators.hpp`.
 6. Available CUDA kernels of the operators are placed in respectives `x_operators.cuh` files. For any CUDA version operator, there is a called `bridge_x` method, managing the exchange between the architectures.

## Licence
The entire code is released under the GPLv3.

## Author

Raphael dos Santos Evangelista and Esteban Walter Gonzalez Clua (Advisor).