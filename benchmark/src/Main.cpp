#include "Output.h"
#include "Parameters.h"
#include "Parser.h"
#include "PerformanceTracker.h"
#include "Simulation.h"

#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

using namespace pica;
using namespace std;


int realMain(int argc, char* argv[]);
long getNumParticleUpdates(const Parameters& parameters);

int main(int argc, char* argv[])
{
    try {
        return realMain(argc, argv);
    }
    catch (std::exception& e) {
        cout << "ERROR: Unhandled exception with message '" << e.what()
            << "', benchmark terminated\n";
        return EXIT_FAILURE;
    }
    catch (...) {
        cout << "ERROR: Unhandled exception, benchmark terminated\n";
        return EXIT_FAILURE;
    }
}

int realMain(int argc, char* argv[])
{
    printHeader();
    Parameters parameters = readParameters(argc, argv);
    printParameters(parameters);
    PerformanceTracker tracker;
    runSimulation(parameters, tracker);
    printPerformance(tracker, getNumParticleUpdates(parameters));
    return 0;
}

long getNumParticleUpdates(const Parameters& parameters)
{
    long numCells = 1;
    for (int d = 0; d < parameters.dimension; d++)
        numCells *= parameters.numCells[d];
    return numCells * parameters.particlesPerCell * parameters.numIterations;
}
