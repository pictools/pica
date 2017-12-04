#include "Output.h"

#include "Parameters.h"
#include "Parser.h"
#include "PerformanceTracker.h"

#include "pica/utility/Utility.h"
#include "pica/threading/OpenMPHelper.h"

#include <iostream>
#include <map>
#include <string>

using namespace pica;
using namespace std;


map<PerformanceTracker::Stage, string> getStageNames();

string getPrefix()
{
    return "   ";
}

void printHeader()
{
    cout << "pica benchmark\n";
}

void printParameters(const Parameters& parameters)
{
    string prefix = getPrefix();
    cout << "\nParameters:\n";
    if (useOpenMP())
        cout << prefix << toString(parameters.numThreads) + " OpenMP threads\n";
    else
        cout << prefix << "OpenMP is disabled\n";
    cout << prefix << "Dimension: " << parameters.dimension << "\n";
    cout << prefix << "Grid size: ";
    for (int d = 0; d < parameters.dimension - 1; d++)
        cout << parameters.numCells[d] << "x";
    cout << parameters.numCells[parameters.dimension - 1] << "\n";
    cout << prefix << "Time iterations: " << parameters.numIterations << "\n";
    cout << prefix << "Particles per cell: " << parameters.particlesPerCell << "\n";
    cout << prefix << "Particles temperature: " << parameters.temperature << "\n";
    cout << prefix << "Particle representation: " << toString(parameters.particleRepresentation) << "\n";
    cout << prefix << "Ensemble representation: " << toString(parameters.ensembleRepresentation);
    if (parameters.ensembleRepresentation == EnsembleRepresentation_Ordered)
        cout << ", sorting period = " << parameters.sortingPeriod;
    else if (parameters.ensembleRepresentation == EnsembleRepresentation_Supercells)
    {
        cout << ", supercell size = " << parameters.numCellsPerSupercell.x << "x" <<
            parameters.numCellsPerSupercell.y << "x" << parameters.numCellsPerSupercell.z;
        cout << ", preloading ";
        if (parameters.enablePreloading)
            cout << "enabled";
        else
            cout << "disabled";
    }
    cout << "\n";
    cout << prefix << "Tile size: " << toString(parameters.tileSize) << "\n";
}

void printPerformance(const PerformanceTracker& tracker, long numParticleUpdates)
{
    string prefix = getPrefix();
    cout << "\nPerformance results:\n";
    map<PerformanceTracker::Stage, string> stageNames = getStageNames();
    PerformanceTracker::StageTime stageTime = tracker.getStageTime();
    double overallTime = 0.0;
    for (PerformanceTracker::StageTime::iterator i = stageTime.begin(); i != stageTime.end(); i++) {
        cout << prefix << stageNames[i->first] << ": " << i->second << " sec.\n";
        overallTime += i->second;
    }
    double nsPerParticleUpdate = overallTime * 1e9 / numParticleUpdates;
    cout << prefix << nsPerParticleUpdate << " ns per particle update\n";
}

map<PerformanceTracker::Stage, string> getStageNames()
{
    map<PerformanceTracker::Stage, string> names;
    names[PerformanceTracker::Stage_FieldSolver] = "Field solver";
    names[PerformanceTracker::Stage_ParticleLoop] = "Particle loop";
    return names;
}
