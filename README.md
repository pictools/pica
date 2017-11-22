# pica

[![Language](https://img.shields.io/badge/language-C%2B%2B-orange.svg)](https://isocpp.org/)
[![Platforms](https://img.shields.io/badge/platform-linux%20%7C%20windows-blue.svg)](../master/)
[![License](https://img.shields.io/badge/license-MIT-lightgrey.svg)](../master/LICENSE)

## About

*pica* is a C++ kernel library for particle-in-cell plasma simulation.

The library provides basic routines and data structures for particle-in-cell plasma simulation. It currently includes Yee grid, first and second-order particle form factors, direct and charge-conserving current deposition, and cell-based particle layout and processing. The main routines are implemented for 1D, 2D and 3D Cartesian coordinates, are optimized and OpenMP-parallelized for multicore CPUs, including Xeon Phi. 

The library is based on the [core of the PICADOR particle-in-cell code](http://www.sciencedirect.com/science/article/pii/S0010465516300194). 

## Licence

*pica* is licensed under the MIT licence. For a detailed description, please refer to [LICENSE](../master/LICENSE).

## How to use

*pica* is currently in the early alpha version and currently contains only implementations for shared memory systems. It is being extended with more data structures and routines, communications and load balancing for distributed memory. We plan to release a beta version along with a skeleton code based on the library on December 2017.
