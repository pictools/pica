#include "Parser.h"

#include "pica/threading/OpenMPHelper.h"
#include "pica/utility/Utility.h"

#include "cmdline.h"

#include <stdexcept>
#include <string>

using namespace pica;
using namespace std;


Parameters readParameters(int argc, char* argv[])
{
    cmdline::parser parser;
    if (useOpenMP())
        parser.add<int>("nthreads", 'n', "number of OpenMP threads", false, getNumThreads());
    parser.add<int>("dimension", 'd', "dimension", false, 3);
    parser.add<int>("grid", 'g', "grid size per each dimension", false, 50);
    parser.add<int>("ppc", 'p', "particles per cell", false, 100);
    parser.add<double>("temperature", 't', "initial temperature", false, 0.0);
    parser.add<int>("iterations", 'i', "number of time iterations", false, 100);
    parser.add<string>("representation_particles", 'r', "representation of particles in arrays: " +
        toString(ParticleRepresentation_SoA) + " or " + toString(ParticleRepresentation_AoS), false,
        toString(ParticleRepresentation_SoA), cmdline::oneof<string>(toString(ParticleRepresentation_SoA),
        toString(ParticleRepresentation_AoS)));
    parser.add<string>("ensemble_storage", 'e', "storage of particles in ensemble: " +
        toString(EnsembleRepresentation_Unordered) + " or " + toString(EnsembleRepresentation_Ordered) +
        " or " + toString(EnsembleRepresentation_Supercells), false, toString(EnsembleRepresentation_Ordered),
        cmdline::oneof<string>(toString(EnsembleRepresentation_Unordered),
        toString(EnsembleRepresentation_Ordered), toString(EnsembleRepresentation_Supercells)));
    parser.parse_check(argc, argv);

    // For now just hardcode
    Parameters parameters;
    parameters.dimension = parser.get<int>("dimension");
    int numCellsEachDimension = parser.get<int>("grid");
    parameters.numCells = Vector3<int>(numCellsEachDimension, numCellsEachDimension, numCellsEachDimension);
    parameters.numIterations = parser.get<int>("iterations");
    if (useOpenMP())
        parameters.numThreads = parser.get<int>("nthreads");
    else
        parameters.numThreads = 1;
    parameters.particlesPerCell = parser.get<int>("ppc");
    parameters.temperature = parser.get<double>("temperature");
    string representationParticles = parser.get<string>("representation_particles");
    if (representationParticles == toString(ParticleRepresentation_SoA))
        parameters.particleRepresentation = ParticleRepresentation_SoA;
    else if (representationParticles == toString(ParticleRepresentation_AoS))
        parameters.particleRepresentation = ParticleRepresentation_AoS;
    else
        throw std::invalid_argument("wrong value of representation particles");
    string ensembleStorage = parser.get<string>("ensemble_storage");
    if (ensembleStorage == toString(EnsembleRepresentation_Unordered))
        parameters.ensembleRepresentation = EnsembleRepresentation_Unordered;
    else if (ensembleStorage == toString(EnsembleRepresentation_Ordered))
        parameters.ensembleRepresentation = EnsembleRepresentation_Ordered;
    else if (ensembleStorage == toString(EnsembleRepresentation_Supercells))
        parameters.ensembleRepresentation = EnsembleRepresentation_Supercells;
    else
        throw std::invalid_argument("wrong value of ensemble representation");
    return parameters;
}
