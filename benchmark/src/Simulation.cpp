#include "Simulation.h"

#include "GridHelper.h"
#include "Parameters.h"
#include "PerformanceTracker.h"

#include "pica/fieldSolver/YeeSolver.h"
#include "pica/grid/Grid.h"
#include "pica/grid/YeeGrid.h"
#include "pica/math/Constants.h"
#include "pica/math/Vectors.h"
#include "pica/threading/OpenMPHelper.h"
#include "pica/utility/Utility.h"

#include <memory>
#include <stdexcept>

using namespace pica;


template<class Grid, class FieldSolver>
void updateField(Grid& grid, FieldSolver& fieldSolver, double dt, PerformanceTracker& tracker)
{
    tracker.start(PerformanceTracker::Stage_FieldSolver);
    fieldSolver.updateB(grid, dt / 2.0);
    fieldSolver.updateE(grid, dt);
    FieldBoundaryConditions<Grid>().apply(grid);
    fieldSolver.updateB(grid, dt / 2.0);
    tracker.finish(PerformanceTracker::Stage_FieldSolver);
}


template<class Grid, class FieldSolver>
void runIteration(Grid& grid, FieldSolver& fieldSolver, double dt, PerformanceTracker& tracker)
{
    updateField(grid, fieldSolver, dt, tracker);
}


template<class Position, class Real>
Real getTimeStep(Position step)
{
    Real sumInvSquares = 0;
    for (int d = 0; d < VectorDimensionHelper<Position>::dimension; d++)
        sumInvSquares += static_cast<Real>(1.0) / (step[d] * step[d]);
    return 0.5 / (Constants<Real>::c() * sqrt(sumInvSquares));
}


template<Dimension dimension>
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
    YeeSolver fieldSolver;
    Real timeStep = getTimeStep<PositionType, Real>(step);
    omp_set_num_threads(parameters.numThreads);
    for (int i = 0; i < parameters.numIterations; i++)
        runIteration(*grid, fieldSolver, timeStep, tracker);
}


void runSimulation(const Parameters& parameters, PerformanceTracker& tracker)
{
    switch (parameters.dimension) {
        case 1:
            runSimulation<One>(parameters, tracker);
            break;
        case 2:
            runSimulation<Two>(parameters, tracker);
            break;
        case 3:
            runSimulation<Three>(parameters, tracker);
            break;
        default:
            throw std::invalid_argument("wrong value of dimension: " + toString(parameters.dimension));
    }
}
