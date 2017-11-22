#ifndef PICA_BENCHMARK_PARAMETERS_H
#define PICA_BENCHMARK_PARAMETERS_H


#include "pica/math/Vectors.h"
#include "pica/particles/Ensemble.h"
#include "pica/particles/ParticleArray.h"


// All parameters are given as in 3D,
// in lower dimensions some components are not used
struct Parameters {
    int dimension;
    pica::Vector3<int> numCells;
    int numIterations;
    int numThreads;
    int particlesPerCell;
    double temperature;
    pica::ParticleRepresentation particleRepresentation;
    pica::EnsembleRepresentation ensembleRepresentation;
    int tileSize;
};


#endif
