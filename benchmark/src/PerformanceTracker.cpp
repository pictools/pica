#include "PerformanceTracker.h"

#include "Timer.h"


PerformanceTracker::PerformanceTracker()
{
    for (int i = 0; i < numStages; i++)
        stageTime[Stage(i)] = 0.0;
}

PerformanceTracker::~PerformanceTracker()
{
    // This is to use std::auto_ptr<Stopwatch> in .h with only forward declaration
}

void PerformanceTracker::start(Stage stage)
{
    stageTimer->reset();
    stageTimer->start();
}

void PerformanceTracker::finish(Stage stage)
{
    stageTimer->stop();
    stageTime[stage] += stageTimer->getElapsed();
}

PerformanceTracker::StageTime PerformanceTracker::getStageTime() const
{
    return stageTime;
}
