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
    parser.add<int>("dimension", 0, "dimension", false, 3);
    parser.add<int>("ncellsx", 0, "number of grid cells in x", false, 40);
    parser.add<int>("ncellsy", 0, "number of grid cells in y, only used for dimension >= 2", false, 40);
    parser.add<int>("ncellsz", 0, "number of grid cells in z, only used for dimension = 3", false, 40);
    parser.add<int>("nppc", 0, "number of particles per cell", false, 50);
    parser.add<double>("temperature", 0, "initial temperature", false, 0.0);
    parser.add<int>("niterations", 0, "number of time iterations", false, 100);
    parser.add<string>("layout", 0, "layout of particles in arrays: " +
        toString(ParticleRepresentation_SoA) + " or " + toString(ParticleRepresentation_AoS), false,
        toString(ParticleRepresentation_SoA), cmdline::oneof<string>(toString(ParticleRepresentation_SoA),
        toString(ParticleRepresentation_AoS)));
    parser.add<string>("ordering", 0, "ordering of particles in ensemble: " +
        toString(EnsembleRepresentation_Unordered) + " or " + toString(EnsembleRepresentation_Ordered) +
        " or " + toString(EnsembleRepresentation_Supercells), false, toString(EnsembleRepresentation_Ordered),
        cmdline::oneof<string>(toString(EnsembleRepresentation_Unordered),
        toString(EnsembleRepresentation_Ordered), toString(EnsembleRepresentation_Supercells)));
    parser.add<int>("sortingperiod", 0, "period of time steps to perform ordering, only used for ordered representation",
        false, 100);
    parser.add<int>("ncellssupercellx", 0, "number of cells per supercells in x, only used for supercells representation",
        false, 2);
    parser.add<int>("ncellssupercelly", 0, "number of cells per supercells in y, only used for supercells representation and dimension >= 2",
        false, 2);
    parser.add<int>("ncellssupercellz", 0, "number of cells per supercells in z, only used for supercells representation and dimension = 3",
        false, 2);
    parser.add("preloading", 0, "enable preloading, only used for supercells ordering");
    parser.add<int>("tilesize", 0, "size of tile for particle processing", false, 8);
    if (useOpenMP())
        parser.add<int>("nthreads", 0, "number of OpenMP threads, default value is based on system settings", false, getNumThreads());
    parser.parse_check(argc, argv);

    // For now just hardcode
    Parameters parameters;
    parameters.dimension = parser.get<int>("dimension");
    parameters.numCells.x = parser.get<int>("ncellsx");
    parameters.numCells.y = parser.get<int>("ncellsy");
    parameters.numCells.z = parser.get<int>("ncellsz");
    parameters.particlesPerCell = parser.get<int>("nppc");
    parameters.temperature = parser.get<double>("temperature");
    parameters.numIterations = parser.get<int>("niterations");
    string representationParticles = parser.get<string>("layout");
    if (representationParticles == toString(ParticleRepresentation_SoA))
        parameters.particleRepresentation = ParticleRepresentation_SoA;
    else if (representationParticles == toString(ParticleRepresentation_AoS))
        parameters.particleRepresentation = ParticleRepresentation_AoS;
    else
        throw std::invalid_argument("wrong value of particle layout");
    string ensembleStorage = parser.get<string>("ordering");
    if (ensembleStorage == toString(EnsembleRepresentation_Unordered))
        parameters.ensembleRepresentation = EnsembleRepresentation_Unordered;
    else if (ensembleStorage == toString(EnsembleRepresentation_Ordered))
        parameters.ensembleRepresentation = EnsembleRepresentation_Ordered;
    else if (ensembleStorage == toString(EnsembleRepresentation_Supercells))
        parameters.ensembleRepresentation = EnsembleRepresentation_Supercells;
    else
        throw std::invalid_argument("wrong value of particle ordering");
    parameters.sortingPeriod = parser.get<int>("sortingperiod");
    parameters.numCellsPerSupercell.x = parser.get<int>("ncellssupercellx");
    parameters.numCellsPerSupercell.y = parser.get<int>("ncellssupercelly");
    parameters.numCellsPerSupercell.z = parser.get<int>("ncellssupercellz");
    parameters.enablePreloading = parser.exist("preloading");
    parameters.tileSize = parser.get<int>("tilesize");
    if (useOpenMP())
        parameters.numThreads = parser.get<int>("nthreads");
    else
        parameters.numThreads = 1;
    return parameters;
}
