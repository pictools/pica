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

void printHeader()
{
    string message = "pica benchmark, ";
    if (useOpenMP())
        message += toString(getNumThreads()) + " OpenMP threads.";
    else
        message += "OpenMP disabled.";
    cout << message << "\n";
}

void printParameters(const Parameters& parameters)
{
    cout << "\nParameters:\n";
    cout << "Dimension = " << parameters.dimension << "\n";
    cout << "Grid size = " << toString(parameters.numCells) << "\n";
    cout << parameters.numIterations << " time steps\n";
}

void printPerformance(const PerformanceTracker& tracker)
{
    cout << "\nPerformance results:\n";
    map<PerformanceTracker::Stage, string> stageNames = getStageNames();
    PerformanceTracker::StageTime stageTime = tracker.getStageTime();
    for (PerformanceTracker::StageTime::iterator i = stageTime.begin(); i != stageTime.end(); i++)
        cout << stageNames[i->first] << ": " << i->second << " sec.\n";
}

map<PerformanceTracker::Stage, string> getStageNames()
{
    map<PerformanceTracker::Stage, string> names;
    names[PerformanceTracker::Stage_CurrentDeposition] = "Current deposition";
    names[PerformanceTracker::Stage_FieldSolver] = "Field solver";
    names[PerformanceTracker::Stage_ParticleLoop] = "Particle loop";
    return names;
}
