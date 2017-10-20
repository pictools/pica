#ifndef PICA_BENCHMARK_PERFORMANCETRACKER_H
#define PICA_BENCHMARK_PERFORMANCETRACKER_H


#include <map>
#include <memory>
#include <vector>


struct Stopwatch;

class PerformanceTracker {
public:
    PerformanceTracker();
    ~PerformanceTracker();

    enum Stage { Stage_CurrentDeposition, Stage_FieldSolver,
        Stage_ParticleLoop, numStages };
    void start(Stage stage);
    void finish(Stage stage);

    typedef std::map<Stage, double> StageTime;
    StageTime getStageTime() const;

private:
    StageTime stageTime;
    std::auto_ptr<Stopwatch> stageTimer;
};


#endif
