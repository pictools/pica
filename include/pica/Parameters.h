#ifndef PICA_PARAMETERS_H
#define PICA_PARAMETERS_H


#include "pica/math/Vectors.h"

#include <iostream>


namespace pica {

class Parameters
{
public:

    FP3 localMin, localMax; // local domain limits
    FP3 globalMin, globalMax; // global simulation area limits
    Int3 localMinIndex, localMaxIndex; // minimum and maximum value of the grid index available in the current domain
    Int3 localGridSize; // number of cells in local domain
    Int3 globalGridSize; // global number of cells
    Int3 gridSize; // MDK alias for globalGridSize
    FP3 step; // coordinate steps
    FP timeStep; // time step
    int numIterations; // Number of interations
    int linearRank; // current MPI process linear identificator
    Int3 rank; // current MPI process 3D identificator
    Int3 numProcesses; // number of MPI processes
    int dimensionality; // number of dimensions in the simulation (1, 2 or 3)

    Parameters():
        timeStep(0),
        numIterations(0),
        linearRank(0),
        numProcesses(1, 1, 1),
        dimensionality(3)
    {}

};

} // namespace pica


#endif
