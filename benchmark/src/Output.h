#ifndef PICA_BENCHMARK_OUTPUT_H
#define PICA_BENCHMARK_OUTPUT_H


struct Parameters;
class PerformanceTracker;

void printHeader();
void printParameters(const Parameters& parameters);
void printPerformance(const PerformanceTracker& tracker, long numParticleUpdates);


#endif
