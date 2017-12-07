#include "Simulation.h"

#include "EnsembleHelper.h"
#include "GridHelper.h"
#include "Parameters.h"
#include "ParticleProcessing.h"
#include "PerformanceTracker.h"

#include "pica/fieldSolver/YeeSolver.h"
#include "pica/grid/Grid.h"
#include "pica/grid/YeeGrid.h"
#include "pica/math/Constants.h"
#include "pica/math/Vectors.h"
#include "pica/particles/Particle.h"
#include "pica/particles/ParticleTraits.h"
#include "pica/threading/OpenMPHelper.h"
#include "pica/utility/Utility.h"

#include <cstdlib>
#include <limits>
#include <memory>
#include <stdexcept>

using namespace pica;


template<class Grid, class FieldSolver>
void updateField(Grid& grid, FieldSolver& fieldSolver, double dt, PerformanceTracker& tracker)
{
    fieldSolver.updateB(grid, dt / 2.0);
    fieldSolver.updateE(grid, dt);
    FieldBoundaryConditions<Grid>().apply(grid);
    fieldSolver.updateB(grid, dt / 2.0);
}


template<class Ensemble, class Grid, class ParticleProcessing, class FieldSolver>
void runIteration(Ensemble& ensemble, Grid& grid, ParticleProcessing& particleProcessing,
    FieldSolver& fieldSolver, double dt, PerformanceTracker& tracker)
{
    tracker.start(PerformanceTracker::Stage_FieldSolver);
    updateField(grid, fieldSolver, dt, tracker);
    tracker.finish(PerformanceTracker::Stage_FieldSolver);

    tracker.start(PerformanceTracker::Stage_ParticleLoop);
    particleProcessing.process(ensemble, grid, dt);
    tracker.finish(PerformanceTracker::Stage_ParticleLoop);
}


template<class Position, class Real>
Real getTimeStep(Position step)
{
    Real sumInvSquares = 0;
    for (int d = 0; d < VectorDimensionHelper<Position>::dimension; d++)
        sumInvSquares += static_cast<Real>(1.0) / (step[d] * step[d]);
    return 0.5 / (Constants<Real>::c() * sqrt(sumInvSquares));
}


double getNormal()
{
    double result;
    do
    {
        double u1 = (double)rand() / (double)RAND_MAX;
        double u2 = (double)rand() / (double)RAND_MAX;
        result = std::sqrt(-2 * std::log(u1)) * std::cos(2 * constants::pi * u2);
    } while (!(result <= std::numeric_limits<double>::max()
        && result >= -std::numeric_limits<double>::max()));
    return result;
}


template<class Ensemble>
void createParticles(const Parameters& parameters, Ensemble& ensemble)
{
    srand(0);
    typedef typename Ensemble::PositionType PositionType;
    typedef typename ParticleTraits<typename Ensemble::Particle>::MomentumType MomentumType;
    int numCells = 1;
    for (int d = 0; d < VectorDimensionHelper<PositionType>::dimension; d++)
        numCells *= parameters.numCells[d];
    int numParticles = numCells * parameters.particlesPerCell;
    PositionType minPosition = ensemble.getMinPosition();
    PositionType maxPosition = ensemble.getMaxPosition();
    for (int i = 0; i < numParticles; i++) {
        typename Ensemble::Particle particle;
        particle.setMass(constants::electronMass);
        particle.setCharge(constants::electronCharge);
        particle.setFactor(1.0);
        PositionType position;
        for (int d = 0; d < VectorDimensionHelper<PositionType>::dimension; d++)
            position[d] = minPosition[d] + (maxPosition[d] - minPosition[d]) *
                (double)rand() / (double)RAND_MAX;
        particle.setPosition(position);
        // The standard deviation is sqrt(1/(2*alpha)), where alpha is
        // 3/2 * ((T/mc^2 + 1)^2 - 1)^(-1)
        double alpha = parameters.temperature / particle.getMass() / constants::c / constants::c + 1;
        alpha = 1.5 / (alpha * alpha - 1);
        double sigma = sqrt(0.5 / alpha) * particle.getMass() * constants::c;
        // Initial particle momentum is combination of given initial
        // momentum based on coords and random term in N(0, sigma)
        MomentumType momentum;
        for (int d = 0; d < VectorDimensionHelper<MomentumType>::dimension; d++)
            momentum[d] = getNormal() * sigma;
        particle.setMomentum(momentum);
        ensemble.add(particle);
    }
}


template<Dimension dimension, class Particle, class ParticleArray, class Ensemble>
void runSimulation(const Parameters& parameters, PerformanceTracker& tracker)
{
    typedef YeeGrid<dimension> Grid;
    typedef typename Grid::PositionType PositionType;
    typedef typename Grid::IndexType IndexType;
    typedef typename Grid::ValueType Real;
    PositionType minPosition;
    PositionType maxPosition = OnesHelper<dimension, Real>::get();
    IndexType numCells;
    for (int d = 0; d < parameters.dimension; d++)
        numCells[d] = parameters.numCells[d];
    PositionType step = (maxPosition - minPosition) / PositionType(numCells);
    std::auto_ptr<Grid> grid = createGrid<Grid>(minPosition, maxPosition, numCells);
    std::auto_ptr<Ensemble> ensemble = createEnsemble<Ensemble>(minPosition, maxPosition,
        numCells, parameters.numCellsPerSupercell);
    createParticles(parameters, *ensemble);
    ParticleProcessing<Ensemble, Grid> particleProcessing(parameters, *ensemble, *grid);
    YeeSolver fieldSolver;
    Real timeStep = getTimeStep<PositionType, Real>(step);
    omp_set_num_threads(parameters.numThreads);
    for (int i = 0; i < parameters.numIterations; i++)
        runIteration(*ensemble, *grid, particleProcessing, fieldSolver, timeStep, tracker);
}

template<Dimension dimension, class Particle, class ParticleArray>
void runSimulation(const Parameters& parameters, PerformanceTracker& tracker)
{
    switch (parameters.ensembleRepresentation) {
        case EnsembleRepresentation_Unordered:
            runSimulation<dimension, Particle, ParticleArray,
                typename Ensemble<ParticleArray, EnsembleRepresentation_Unordered>::Type>(parameters, tracker);
            break;
        case EnsembleRepresentation_Ordered:
            runSimulation<dimension, Particle, ParticleArray,
                typename Ensemble<ParticleArray, EnsembleRepresentation_Ordered>::Type>(parameters, tracker);
            break;
        case EnsembleRepresentation_Supercells:
            runSimulation<dimension, Particle, ParticleArray,
                typename Ensemble<ParticleArray, EnsembleRepresentation_Supercells>::Type>(parameters, tracker);
            break;
        default:
            throw std::invalid_argument("wrong value of ensemble representation");
    }
}


template<Dimension dimension, class Particle>
void runSimulation(const Parameters& parameters, PerformanceTracker& tracker)
{
    switch (parameters.particleRepresentation) {
        case ParticleRepresentation_AoS:
            runSimulation<dimension, Particle, typename ParticleArray<Particle, ParticleRepresentation_AoS>::Type>(parameters, tracker);
            break;
        case ParticleRepresentation_SoA:
            runSimulation<dimension, Particle, typename ParticleArray<Particle, ParticleRepresentation_SoA>::Type>(parameters, tracker);
            break;
        default:
            throw std::invalid_argument("wrong value of particle representation");
    }
}


void runSimulation(const Parameters& parameters, PerformanceTracker& tracker)
{
    switch (parameters.dimension) {
        /*case 1:
            runSimulation<One, Particle1d>(parameters, tracker);
            break;
        case 2:
            runSimulation<Two, Particle2d>(parameters, tracker);
            break;*/
        case 3:
            runSimulation<Three, Particle3d>(parameters, tracker);
            break;
        default:
            throw std::invalid_argument("wrong value of dimension: " + toString(parameters.dimension));
    }
}
