#ifndef PICA_BENCHMARK_PARTICLEPROCESSING_H
#define PICA_BENCHMARK_PARTICLEPROCESSING_H


#include "Parameters.h"

#include "pica/currentDeposition/CurrentDepositor.h"
#include "pica/fieldInterpolation/FieldInterpolator.h"
#include "pica/math/Vectors.h"
#include "pica/particles/Ensemble.h"
#include "pica/particles/ParticleTraits.h"
#include "pica/particlePush/BorisPusher.h"
#include "pica/threading/OpenMPHelper.h"

#include <algorithm>
#include <memory>


template<class Ensemble, class Grid>
struct ParticleProcessing {
    ParticleProcessing(const Parameters&);
    void process(Ensemble& ensemble, Grid& grid, double dt);
};


namespace internal {

template<class Ensemble, class ParticleArray, class Grid, int tileSize>
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
        const int numParticles = endIdx - beginIdx;
        const int numTiles = (numParticles + tileSize - 1) / tileSize;
        for (int tileIdx = 0; tileIdx < numTiles; tileIdx++) {
            const int tileBeginIdx = beginIdx + tileIdx * tileSize;
            const int tileEndIdx = std::min(tileBeginIdx + tileSize, endIdx);
            processTile(particles, tileBeginIdx, tileEndIdx, grid, dt);
        }
    }

private:
    PositionType minPosition, maxPosition;

    void processTile(ParticleArray& particles, int beginIdx, int endIdx, Grid& grid, double dt)
    {
        pushParticles(particles, beginIdx, endIdx, grid, dt);
        applyBoundaryConditions(particles, beginIdx, endIdx, grid, dt);
        depositCurrents(particles, beginIdx, endIdx, grid, dt);
    }

    void pushParticles(ParticleArray& particles, int beginIdx, int endIdx, Grid& grid, double dt)
    {
        MomentumType interpolatedE[tileSize];
        MomentumType interpolatedB[tileSize];
        pica::FieldInterpolatorCIC<Grid> fieldInterpolator(grid);
        for (int particleIdx = beginIdx; particleIdx < endIdx; particleIdx++)
            fieldInterpolator.get(particles[particleIdx].getPosition(),
                interpolatedE[particleIdx - beginIdx], interpolatedB[particleIdx - beginIdx]);

        pica::BorisPusher pusher;
        for (int particleIdx = beginIdx; particleIdx < endIdx; particleIdx++)
            pusher.push<ParticleRef, MomentumType, PositionType, double>(particles[particleIdx],
                interpolatedE[particleIdx - beginIdx], interpolatedB[particleIdx - beginIdx], dt);
    }

    void applyBoundaryConditions(ParticleArray& particles, int beginIdx, int endIdx, Grid& grid, double dt)
    {
        // Reflecting boundary conditions
        for (int particleIdx = beginIdx; particleIdx < endIdx; particleIdx++) {
            PositionType position = particles[particleIdx].getPosition();
            MomentumType momentum = particles[particleIdx].getMomentum();
            for (int d = 0; d < VectorDimensionHelper<PositionType>::dimension; d++)
                if (position[d] < minPosition[d]) {
                    position[d] = 2 * minPosition[d] - position[d];
                    momentum[d] = -momentum[d];
                }
                else
                    if (position[d] > maxPosition[d]) {
                        position[d] = 2 * maxPosition[d] - position[d];
                        momentum[d] = -momentum[d];
                    }
            particles[particleIdx].setPosition(position);
            particles[particleIdx].setMomentum(momentum);
        }
    }

    void depositCurrents(ParticleArray& particles, int beginIdx, int endIdx, Grid& grid, double dt)
    {
        const double halfDt = 0.5 * dt;
        pica::CurrentDepositorCIC<Grid> currentDepositor(grid);
        for (int particleIdx = beginIdx; particleIdx < endIdx; particleIdx++) {
            PositionType position = particles[particleIdx].getPosition();
            for (int d = 0; d < VectorDimensionHelper<PositionType>::dimension; d++)
                position[d] -= particles[particleIdx].getVelocity()[d] * halfDt;
            MomentumType current = particles[particleIdx].getVelocity() * particles[particleIdx].getCharge() * (double)particles[particleIdx].getFactor();
            currentDepositor.deposit(position, current);
        }
    }

};

} // namespace internal


template<class ParticleArray, class Grid>
struct ParticleProcessing<pica::EnsembleUnordered<ParticleArray>, Grid> {
    typedef pica::EnsembleUnordered<ParticleArray> Ensemble;

    ParticleProcessing(const Parameters& parameters) :
        pImpl(createImplementation(parameters))
    {
    }

    void process(Ensemble& ensemble, Grid& grid, double dt)
    {
        pImpl->process(ensemble, grid, dt);
    }

private:

    class ImplementationBase {
    public:
        virtual ~ImplementationBase() {}
        virtual void process(Ensemble& ensemble, Grid& grid, double dt) = 0;
    };

    template<int tileSize>
    class Implementation : public ImplementationBase {
    public:

        virtual void process(Ensemble& ensemble, Grid& grid, double dt)
        {
            zeroizeCurrents(grid);
            processParticles(ensemble, grid, dt);
            finalizeCurrents(grid);
        }

    private:

        void zeroizeCurrents(Grid& grid)
        {
        }

        void processParticles(Ensemble& ensemble, Grid& grid, double dt)
        {
            typedef ::internal::ParticleArrayProcessing<Ensemble, ParticleArray, Grid, tileSize> ParticleArrayProcessing;
            ParticleArrayProcessing particleArrayProcessing(ensemble);
            const int numParticles = ensemble.size();
            const int numThreads = getNumThreads();
            const int particlesPerThread = (numParticles + numThreads - 1) / numThreads;
            #pragma omp parallel for
            for (int idx = 0; idx < numThreads; idx++) {
                const int beginIdx = idx * particlesPerThread;
                const int endIdx = std::min(beginIdx + particlesPerThread, numParticles);
                particleArrayProcessing.process(ensemble.getParticles(), beginIdx, endIdx, grid, dt);
            }
        }

        void finalizeCurrents(Grid& grid)
        {
        }

    };

    std::auto_ptr<ImplementationBase> pImpl;

    ImplementationBase* createImplementation(const Parameters& parameters)
    {
        const int tileSize = parameters.tileSize;
        if (tileSize <= 8)
            return new Implementation<8>;
        else if (tileSize <= 16)
            return new Implementation<16>;
        else if (tileSize <= 24)
            return new Implementation<24>;
        else if (tileSize <= 32)
            return new Implementation<32>;
        else if (tileSize <= 40)
            return new Implementation<40>;
        else if (tileSize <= 48)
            return new Implementation<48>;
        else if (tileSize <= 56)
            return new Implementation<56>;
        else if (tileSize <= 64)
            return new Implementation<64>;
        else if (tileSize <= 72)
            return new Implementation<72>;
        else if (tileSize <= 80)
            return new Implementation<80>;
        else if (tileSize <= 88)
            return new Implementation<88>;
        else if (tileSize <= 96)
            return new Implementation<96>;
        else if (tileSize <= 104)
            return new Implementation<104>;
        else if (tileSize <= 112)
            return new Implementation<112>;
        else if (tileSize <= 120)
            return new Implementation<120>;
        else
            return new Implementation<128>;
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
