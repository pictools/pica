#include "Output.h"

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
        cout << prefix << toString(getNumThreads()) + " OpenMP threads\n";
    else
        cout << prefix << "OpenMP is disabled\n";
    cout << prefix << "Dimension: " << parameters.dimension << "\n";
    cout << prefix << "Grid size: " << toString(parameters.numCells) << "\n";
    cout << prefix << "Time iterations: " << parameters.numIterations << "\n";
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
