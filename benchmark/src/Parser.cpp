#include "Parser.h"

#include "pica/math/Constants.h"

using namespace pica;


Parameters readParameters(int argc, char* argv[])
{
    // For now just hardcode
    Parameters parameters;
    parameters.dimension = 3;
    parameters.numCells = Vector3<int>(32, 32, 32);
    parameters.numIterations = 100;
    return parameters;
}
