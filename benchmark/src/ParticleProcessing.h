#ifndef PICA_BENCHMARK_PARTICLEPROCESSING_H
#define PICA_BENCHMARK_PARTICLEPROCESSING_H


#include "Parameters.h"

#include "pica/fieldInterpolation/FieldInterpolator.h"
#include "pica/math/Vectors.h"
#include "pica/particles/Ensemble.h"
#include "pica/particles/ParticleTraits.h"
#include "pica/particlePush/BorisPusher.h"
#include "pica/threading/OpenMPHelper.h"

#include <algorithm>


template<class Ensemble, class Grid>
struct ParticleProcessing {
    ParticleProcessing(const Parameters&);
    void process(Ensemble& ensemble, Grid& grid, double dt);
};


namespace internal {

template<class Ensemble, class ParticleArray, class Grid>
class ParticleArrayProcessing {
public:
    typedef typename Ensemble::Particle Particle;
    typedef typename Ensemble::ParticleRef ParticleRef;
    typedef typename pica::ParticleTraits<Particle>::MomentumType MomentumType;
    typedef typename pica::ParticleTraits<Particle>::PositionType PositionType;

    ParticleArrayProcessing(const Ensemble& ensemble) :
        minPosition(ensemble.getMinPosition()),
        maxPosition(ensemble.getMaxPosition())
    {}

    void process(ParticleArray& particles, int beginIdx, int endIdx, Grid& grid, double dt)
    {
        pica::BorisPusher pusher;
        pica::FieldInterpolatorCIC<Grid> fieldInterpolator(grid);
        for (int i = beginIdx; i < endIdx; i++) {
            MomentumType e, b;
            fieldInterpolator.get(particles[i].getPosition(), e, b);
            pusher.push<ParticleRef, MomentumType, PositionType, double>(particles[i], e, b, dt);

            // reflecting boundary conditions
            PositionType position = particles[i].getPosition();
            MomentumType momentum = particles[i].getMomentum();
            for (int d = 0; d < VectorDimensionHelper<PositionType>::dimension; d++)
                if (position[d] < minPosition[d]) {
                    position[d] = 2 * minPosition[d] - position[d];
                    momentum[d] = -momentum[d];
                }
                else if (position[d] > maxPosition[d]) {
                    position[d] = 2 * maxPosition[d] - position[d];
                    momentum[d] = -momentum[d];
                }
            particles[i].setPosition(position);
            particles[i].setMomentum(momentum);
        }
    }

private:
    PositionType minPosition, maxPosition;

};

}


template<class ParticleArray, class Grid>
struct ParticleProcessing<pica::EnsembleUnordered<ParticleArray>, Grid> {
    typedef pica::EnsembleUnordered<ParticleArray> Ensemble;

    ParticleProcessing(const Parameters& ) {}

    void process(Ensemble& ensemble, Grid& grid, double dt)
    {
        internal::ParticleArrayProcessing<Ensemble, ParticleArray, Grid> particleArrayProcessing(ensemble);
        process(particleArrayProcessing, ensemble, grid, dt);
    }

private:

    template<class ParticleArrayProcessing>
    void process(ParticleArrayProcessing& particleArrayProcessing, Ensemble& ensemble, Grid& grid, double dt)
    {
        const int numParticles = ensemble.size();
        const int numTiles = getNumThreads();
        const int tileSize = (numParticles + numTiles - 1) / numTiles;
        #pragma omp parallel for
        for (int tileIdx = 0; tileIdx < numTiles; tileIdx++) {
            const int beginIdx = tileIdx * tileSize;
            const int endIdx = std::min(beginIdx + tileSize, numParticles);
            particleArrayProcessing.process(ensemble.getParticles(), beginIdx, endIdx, grid, dt);
        }
    }

};


template<class ParticleArray, class Grid>
struct ParticleProcessing<pica::EnsembleOrdered<ParticleArray>, Grid> : 
        ParticleProcessing<pica::EnsembleUnordered<ParticleArray>, Grid> {

    ParticleProcessing(const Parameters& parameters) :
        ParticleProcessing<pica::EnsembleUnordered<ParticleArray>, Grid>(parameters),
        iteration(0),
        sortingPeriod(100)
    {}

    void process(pica::EnsembleOrdered<ParticleArray>& ensemble, Grid& grid, double dt)
    {
        if (iteration % sortingPeriod == 0)
            ensemble.reorder();
        iteration++;
        ParticleProcessing<pica::EnsembleUnordered<ParticleArray>, Grid>::process(ensemble, grid, dt);
    }

private:
    int iteration, sortingPeriod;
};


#endif
