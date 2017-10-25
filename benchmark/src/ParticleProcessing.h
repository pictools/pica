#ifndef PICA_BENCHMARK_PARTICLEPROCESSING_H
#define PICA_BENCHMARK_PARTICLEPROCESSING_H


#include "Parameters.h"

#include "pica/fieldInterpolation/FieldInterpolator.h"
#include "pica/math/Vectors.h"
#include "pica/particles/Ensemble.h"
#include "pica/particles/ParticleTraits.h"
#include "pica/particlePush/BorisPusher.h"
#include "pica/threading/OpenMPHelper.h"


template<class Ensemble, class Grid>
struct ParticleProcessing {
    ParticleProcessing(const Parameters&);
    void process(Ensemble& ensemble, Grid& grid, double dt);
};


template<class ParticleArray, class Grid>
struct ParticleProcessing<pica::EnsembleUnordered<ParticleArray>, Grid> {
    typedef pica::EnsembleUnordered<ParticleArray> Ensemble;
    typedef typename Ensemble::Particle Particle;
    typedef typename Ensemble::ParticleRef ParticleRef;
    typedef typename pica::ParticleTraits<Particle>::MomentumType MomentumType;
    typedef typename pica::ParticleTraits<Particle>::PositionType PositionType;

    ParticleProcessing(const Parameters&) {}
    void process(Ensemble& ensemble, Grid& grid, double dt)
    {
        pica::BorisPusher pusher;
        pica::FieldInterpolatorCIC<Grid> fieldInterpolator(grid);
        PositionType minPosition = ensemble.getMinPosition();
        PositionType maxPosition = ensemble.getMaxPosition();
        const int numParticles = ensemble.size();
        #pragma omp parallel for
        for (int i = 0; i < numParticles; i++) {
            MomentumType e, b;
            fieldInterpolator.get(ensemble[i].getPosition(), e, b);
            pusher.push<ParticleRef, MomentumType, PositionType, double>(ensemble[i], e, b, dt);

            // reflecting boundary conditions
            PositionType position = ensemble[i].getPosition();
            MomentumType momentum = ensemble[i].getMomentum();
            for (int d = 0; d < VectorDimensionHelper<PositionType>::dimension; d++)
                if (position[d] < minPosition[d]) {
                    position[d] = 2 * minPosition[d] - position[d];
                    momentum[d] = -momentum[d];
                }
                else if (position[d] > maxPosition[d]) {
                    position[d] = 2 * maxPosition[d] - position[d];
                    momentum[d] = -momentum[d];
                }
            ensemble[i].setPosition(position);
            ensemble[i].setMomentum(momentum);
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
