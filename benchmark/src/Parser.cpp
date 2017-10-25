#include "Parser.h"

#include "pica/threading/OpenMPHelper.h"

using namespace pica;


Parameters readParameters(int argc, char* argv[])
{
    // For now just hardcode
    Parameters parameters;
    parameters.dimension = 3;
    parameters.numCells = Vector3<int>(32, 32, 32);
    parameters.numIterations = 100;
    parameters.numThreads = getNumThreads();
    parameters.particlesPerCell = 30;
    parameters.temperature = 0.0;
    parameters.particleRepresentation = ParticleRepresentation_SoA;
    parameters.ensembleRepresentation = EnsembleRepresentation_Unordered;
    return parameters;
}
