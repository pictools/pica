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

string toString(ParticleRepresentation particleRepresentation)
{
    map<ParticleRepresentation, string> names;
    names[ParticleRepresentation_AoS] = "AoS";
    names[ParticleRepresentation_SoA] = "SoA";
    return names[particleRepresentation];
}

string toString(EnsembleRepresentation ensembleRepresentation)
{
    map<EnsembleRepresentation, string> names;
    names[EnsembleRepresentation_Unordered] = "unordered";
    names[EnsembleRepresentation_Ordered] = "ordered";
    return names[ensembleRepresentation];
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
    cout << prefix << "Ensemble representation: " << toString(parameters.ensembleRepresentation) << "\n";
}

void printPerformance(const PerformanceTracker& tracker)
{
    string prefix = getPrefix();
    cout << "\nPerformance results:\n";
    map<PerformanceTracker::Stage, string> stageNames = getStageNames();
    PerformanceTracker::StageTime stageTime = tracker.getStageTime();
    for (PerformanceTracker::StageTime::iterator i = stageTime.begin(); i != stageTime.end(); i++)
        cout << prefix << stageNames[i->first] << ": " << i->second << " sec.\n";
}

map<PerformanceTracker::Stage, string> getStageNames()
{
    map<PerformanceTracker::Stage, string> names;
    names[PerformanceTracker::Stage_CurrentDeposition] = "Current deposition";
    names[PerformanceTracker::Stage_FieldSolver] = "Field solver";
    names[PerformanceTracker::Stage_ParticleLoop] = "Particle loop";
    return names;
}
