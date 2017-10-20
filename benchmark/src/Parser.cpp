#include "Parser.h"

#include "pica/math/Constants.h"

using namespace pica;


Parameters readParameters(int argc, char* argv[])
{
    // For now just hardcode
    Parameters parameters;
    parameters.dimension = 3;
    parameters.minPosition = Vector3<double>(0.0, 0.0, 0.0);
    parameters.maxPosition = Vector3<double>(1.0, 1.0, 1.0);
    parameters.numCells = Vector3<int>(32, 32, 32);
    parameters.numIterations = 100;
    parameters.dt = 1e-2 / Constants<double>::c();
    return parameters;
}
