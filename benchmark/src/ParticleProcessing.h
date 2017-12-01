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
#include <vector>


template<class Ensemble, class Grid>
struct ParticleProcessing {
    ParticleProcessing(const Parameters&);
    void process(Ensemble& ensemble, Grid& grid, double dt);
};


namespace internal {

// This class implements processing of array of particles
// Processing is tiled according to given tileSize
template<class Ensemble, class ParticleArray, int tileSize>
class ParticleArrayProcessing {
public:
    typedef typename Ensemble::Particle Particle;
    typedef typename pica::ParticleTraits<Particle>::MomentumType MomentumType;
    typedef typename pica::ParticleTraits<Particle>::PositionType PositionType;

    template<class ParticlePusher, class FieldInterpolator, class CurrentDepositor,
        class BoundaryConditions>
    static void process(ParticleArray& particles, int beginIdx, int endIdx,
        ParticlePusher& particlePusher, FieldInterpolator& fieldInterpolator,
        CurrentDepositor& currentDepositor, BoundaryConditions& boundaryConditions, double dt)
    {
        const int numParticles = endIdx - beginIdx;
        const int numTiles = (numParticles + tileSize - 1) / tileSize;
        for (int tileIdx = 0; tileIdx < numTiles; tileIdx++) {
            const int tileBeginIdx = beginIdx + tileIdx * tileSize;
            const int tileEndIdx = std::min(tileBeginIdx + tileSize, endIdx);
            processTile(particles, tileBeginIdx, tileEndIdx, particlePusher, fieldInterpolator,
                currentDepositor, boundaryConditions, dt);
        }
    }

private:

    template<class ParticlePusher, class FieldInterpolator, class CurrentDepositor,
        class BoundaryConditions>
    static void processTile(ParticleArray& particles, int beginIdx, int endIdx,
        ParticlePusher& particlePusher, FieldInterpolator& fieldInterpolator,
        CurrentDepositor& currentDepositor, BoundaryConditions& boundaryConditions, double dt)
    {
        pushParticles(particles, beginIdx, endIdx, particlePusher, fieldInterpolator, dt);
        applyBoundaryConditions(particles, beginIdx, endIdx, boundaryConditions);
        depositCurrents(particles, beginIdx, endIdx, currentDepositor, dt);
    }

    template<class ParticlePusher, class FieldInterpolator>
    static void pushParticles(ParticleArray& particles, int beginIdx, int endIdx,
        ParticlePusher& particlePusher, FieldInterpolator& fieldInterpolator, double dt)
    {
        MomentumType interpolatedE[tileSize];
        MomentumType interpolatedB[tileSize];
        for (int particleIdx = beginIdx; particleIdx < endIdx; particleIdx++)
            fieldInterpolator.get(particles[particleIdx].getPosition(),
                interpolatedE[particleIdx - beginIdx], interpolatedB[particleIdx - beginIdx]);

        for (int particleIdx = beginIdx; particleIdx < endIdx; particleIdx++)
            particlePusher.push(particles[particleIdx],
                interpolatedE[particleIdx - beginIdx], interpolatedB[particleIdx - beginIdx], dt);
    }

    template<class BoundaryConditions>
    static void applyBoundaryConditions(ParticleArray& particles, int beginIdx, int endIdx,
        BoundaryConditions& boundaryConditions)
    {
        for (int particleIdx = beginIdx; particleIdx < endIdx; particleIdx++)
            boundaryConditions.apply(particles[particleIdx]);
    }

    template<class CurrentDepositor>
    static void depositCurrents(ParticleArray& particles, int beginIdx, int endIdx,
        CurrentDepositor& currentDepositor, double dt)
    {
        const double halfDt = 0.5 * dt;
        for (int particleIdx = beginIdx; particleIdx < endIdx; particleIdx++) {
            PositionType position = particles[particleIdx].getPosition();
            for (int d = 0; d < VectorDimensionHelper<PositionType>::dimension; d++)
                position[d] -= particles[particleIdx].getVelocity()[d] * halfDt;
            MomentumType current = particles[particleIdx].getVelocity() * particles[particleIdx].getCharge() * (double)particles[particleIdx].getFactor();
            currentDepositor.deposit(position, current);
        }
    }

};


template<class ParticleArray>
class ReflectingBoundaryConditions {
public:

    typedef typename ParticleArray::Particle Particle;
    typedef typename pica::ParticleTraits<Particle>::MomentumType MomentumType;
    typedef typename pica::ParticleTraits<Particle>::PositionType PositionType;

    ReflectingBoundaryConditions(const PositionType& minPosition, const PositionType& maxPosition) :
        minPosition(minPosition),
        maxPosition(maxPosition)
    {}

    template<class ParticleRef>
    void apply(ParticleRef particle)
    {
        PositionType position = particle.getPosition();
        MomentumType momentum = particle.getMomentum();
        for (int d = 0; d < VectorDimensionHelper<PositionType>::dimension; d++)
            if (position[d] < minPosition[d]) {
                position[d] = 2.0 * minPosition[d] - position[d];
                momentum[d] = -momentum[d];
            }
            else
                if (position[d] > maxPosition[d]) {
                    position[d] = 2.0 * maxPosition[d] - position[d];
                    momentum[d] = -momentum[d];
                }
        particle.setPosition(position);
        particle.setMomentum(momentum);
    }

private:

    PositionType minPosition, maxPosition;

};


} // namespace internal


template<class ParticleArray, class Grid>
struct ParticleProcessing<pica::EnsembleUnordered<ParticleArray>, Grid> {
    typedef pica::EnsembleUnordered<ParticleArray> Ensemble;
    typedef typename Ensemble::Particle Particle;

    ParticleProcessing(const Parameters& parameters, const Ensemble& ensemble, const Grid& grid) :
        pImpl(createImplementation(parameters))
    {
        pImpl->init(parameters, ensemble, grid);
    }

    void process(Ensemble& ensemble, Grid& grid, double dt)
    {
        pImpl->process(ensemble, grid, dt);
    }

private:

    class ImplementationBase {
    public:
        virtual void init(const Parameters& parameters, const Ensemble& ensemble, const Grid& grid) {}
        virtual ~ImplementationBase() {}
        virtual void process(Ensemble& ensemble, Grid& grid, double dt) = 0;
    };

    template<int tileSize>
    class Implementation : public ImplementationBase {
    public:

        virtual void init(const Parameters& parameters, const Ensemble& ensemble, const Grid& grid)
        {
            threadGridCopies.resize(pica::getNumThreads(), grid);
        }

        virtual void process(Ensemble& ensemble, Grid& grid, double dt)
        {
            zeroizeCurrents();
            processParticles(ensemble, grid, dt);
            finalizeCurrents(grid);
        }

    private:

        void zeroizeCurrents()
        {
            #pragma omp parallel for
            for (int threadIdx = 0; threadIdx < threadGridCopies.size(); threadIdx++) {
                Grid& grid = threadGridCopies[threadIdx];
                // The following is for 3D only
                for (int i = 0; i < grid.getSize().x; i++)
                for (int j = 0; j < grid.getSize().y; j++)
                for (int k = 0; k < grid.getSize().z; k++) {
                    grid.jx(i, j, k) = 0.0;
                    grid.jy(i, j, k) = 0.0;
                    grid.jz(i, j, k) = 0.0;
                }
            }
        }

        void processParticles(Ensemble& ensemble, Grid& grid, double dt)
        {
            typedef ::internal::ParticleArrayProcessing<Ensemble, ParticleArray, tileSize> ParticleArrayProcessing;
            typedef ::internal::ReflectingBoundaryConditions<ParticleArray> BoundaryConditions;
            typedef pica::BorisPusher<Particle> ParticlePusher;
            const int numParticles = ensemble.size();
            const int numThreads = getNumThreads();
            const int particlesPerThread = (numParticles + numThreads - 1) / numThreads;
            #pragma omp parallel for
            for (int idx = 0; idx < numThreads; idx++) {
                ParticlePusher particlePusher;
                pica::FieldInterpolatorCIC<Grid> fieldInterpolator(grid);
                pica::CurrentDepositorCIC<Grid> currentDepositor(threadGridCopies[omp_get_thread_num()]);
                BoundaryConditions boundaryConditions(ensemble.getMinPosition(), ensemble.getMaxPosition());
                const int beginIdx = idx * particlesPerThread;
                const int endIdx = std::min(beginIdx + particlesPerThread, numParticles);
                ParticleArrayProcessing::process(ensemble.getParticles(), beginIdx, endIdx,
                    particlePusher, fieldInterpolator, currentDepositor, boundaryConditions, dt);
            }
        }

        void finalizeCurrents(Grid& grid)
        {
            // The following is for 3D only
            double normalization = 1.0 / grid.getStep().volume();
            #pragma omp parallel for collapse(2)
            for (int i = 0; i < grid.getSize().x; i++)
            for (int j = 0; j < grid.getSize().y; j++)
            for (int k = 0; k < grid.getSize().z; k++) {
                grid.jx(i, j, k) = 0.0;
                grid.jy(i, j, k) = 0.0;
                grid.jz(i, j, k) = 0.0;
                for (int threadIdx = 0; threadIdx < threadGridCopies.size(); threadIdx++) {
                    const Grid& currentGrid = threadGridCopies[threadIdx];
                    grid.jx(i, j, k) += currentGrid.jx(i, j, k);
                    grid.jy(i, j, k) += currentGrid.jy(i, j, k);
                    grid.jz(i, j, k) += currentGrid.jz(i, j, k);
                }
                grid.jx(i, j, k) *= normalization;
                grid.jy(i, j, k) *= normalization;
                grid.jz(i, j, k) *= normalization;
            }
        }

        std::vector<Grid> threadGridCopies;

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

    ParticleProcessing(const Parameters& parameters, const Ensemble& ensemble, const Grid& grid) :
        ParticleProcessing<pica::EnsembleUnordered<ParticleArray>, Grid>(parameters, ensemble, grid),
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
