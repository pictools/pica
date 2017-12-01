#ifndef PICA_BENCHMARK_GRIDHELPER_H
#define PICA_BENCHMARK_GRIDHELPER_H


#include "pica/grid/YeeGrid.h"
#include "pica/math/Dimension.h"
#include "pica/math/Vectors.h"
#include "pica/threading/OpenMPHelper.h"

#include <memory>


// Number of ghost cells for each side and dimension
inline int getNumGhostCells()
{
    return 2;
}

template<class Grid>
std::auto_ptr<Grid> createGrid(typename Grid::PositionType minPosition,
    typename Grid::PositionType maxPosition, typename Grid::IndexType numInternalCells)
{
    typedef typename Grid::PositionType PositionType;
    typedef typename Grid::IndexType IndexType;
    PositionType step = (maxPosition - minPosition) / Grid::PositionType(numInternalCells);
    int numGhostCells = getNumGhostCells();
    PositionType origin = minPosition - step * static_cast<typename pica::ScalarType<PositionType>::Type>(numGhostCells);
    IndexType numCells = numInternalCells;
    for (int d = 0; d < pica::VectorDimensionHelper<IndexType>::dimension; d++)
        numCells[d] += 2 * numGhostCells;
    return std::auto_ptr<Grid>(new Grid(origin, step, numCells));
}

template<class Grid>
static typename Grid::IndexType getNumInternalCells(Grid& grid, int numGhostCells)
{
    /// Logic is hardcoded for now to match YeeSolver
    int dimension = pica::VectorDimensionHelper<typename Grid::IndexType>::dimension;
    typename Grid::IndexType numInternalCells = grid.getSize();
    for (int d = 0; d < dimension; d++)
        numInternalCells[d] -= 2 * numGhostCells;
    return numInternalCells;
}

template<class Grid>
struct FieldBoundaryConditions {
    void apply(Grid& grid);
};

template<class Real>
struct FieldBoundaryConditions<pica::YeeGrid<pica::One, Real>> {
    void apply(pica::YeeGrid<pica::One, Real>& grid)
    {
        typedef typename pica::YeeGrid<pica::One, Real>::IndexType IndexType;
        int numGhostCells = getNumGhostCells();
        IndexType numInternalCells = getNumInternalCells(grid, numGhostCells);
        for (int i = 0; i < numGhostCells; i++) {
            grid.ex(i) = grid.ex(i + numInternalCells.x);
            grid.ey(i) = grid.ey(i + numInternalCells.x);
            grid.ez(i) = grid.ez(i + numInternalCells.x);
            grid.ex(numInternalCells.x + numGhostCells + i) = grid.ex(numGhostCells + i);
            grid.ey(numInternalCells.x + numGhostCells + i) = grid.ey(numGhostCells + i);
            grid.ez(numInternalCells.x + numGhostCells + i) = grid.ez(numGhostCells + i);
        }
    }
};

template<class Real>
struct FieldBoundaryConditions<pica::YeeGrid<pica::Two, Real>> {
    void apply(pica::YeeGrid<pica::Two, Real>& grid)
    {
        typedef typename pica::YeeGrid<pica::Two, Real>::IndexType IndexType;
        int numGhostCells = getNumGhostCells();
        IndexType numCells = grid.getSize();
        IndexType numInternalCells = getNumInternalCells(grid, numGhostCells);
        #pragma omp parallel for
        for (int j = 0; j < numCells.y; j++)
            for (int i = 0; i < numGhostCells; i++) {
                grid.ex(i, j) = grid.ex(i + numInternalCells.x, j);
                grid.ey(i, j) = grid.ey(i + numInternalCells.x, j);
                grid.ez(i, j) = grid.ez(i + numInternalCells.x, j);
                grid.ex(numInternalCells.x + numGhostCells + i, j) = grid.ex(numGhostCells + i, j);
                grid.ey(numInternalCells.x + numGhostCells + i, j) = grid.ey(numGhostCells + i, j);
                grid.ez(numInternalCells.x + numGhostCells + i, j) = grid.ez(numGhostCells + i, j);
            }
        #pragma omp parallel for
        for (int i = 0; i < numCells.x; i++)
            for (int j = 0; j < numGhostCells; j++) {
                grid.ex(i, j) = grid.ex(i, j + numInternalCells.y);
                grid.ey(i, j) = grid.ey(i, j + numInternalCells.y);
                grid.ez(i, j) = grid.ez(i, j + numInternalCells.y);
                grid.ex(i, numInternalCells.y + numGhostCells + j) = grid.ex(i, numGhostCells + j);
                grid.ey(i, numInternalCells.y + numGhostCells + j) = grid.ey(i, numGhostCells + j);
                grid.ez(i, numInternalCells.y + numGhostCells + j) = grid.ez(i, numGhostCells + j);
            }
    }
};

template<class Real>
struct FieldBoundaryConditions<pica::YeeGrid<pica::Three, Real>> {
    void apply(pica::YeeGrid<pica::Three, Real>& grid)
    {
        typedef typename pica::YeeGrid<pica::Three, Real>::IndexType IndexType;
        int numGhostCells = getNumGhostCells();
        IndexType numCells = grid.getSize();
        IndexType numInternalCells = getNumInternalCells(grid, numGhostCells);
        #pragma omp parallel for collapse(2)
        for (int j = 0; j < numCells.y; j++)
            for (int k = 0; k < numCells.z; k++)
                for (int i = 0; i < numGhostCells; i++) {
                    grid.ex(i, j, k) = grid.ex(i + numInternalCells.x, j, k);
                    grid.ey(i, j, k) = grid.ey(i + numInternalCells.x, j, k);
                    grid.ez(i, j, k) = grid.ez(i + numInternalCells.x, j, k);
                    grid.ex(numInternalCells.x + numGhostCells + i, j, k) = grid.ex(numGhostCells + i, j, k);
                    grid.ey(numInternalCells.x + numGhostCells + i, j, k) = grid.ey(numGhostCells + i, j, k);
                    grid.ez(numInternalCells.x + numGhostCells + i, j, k) = grid.ez(numGhostCells + i, j, k);
                }
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < numCells.x; i++)
            for (int k = 0; k < numCells.z; k++)
                for (int j = 0; j < numGhostCells; j++) {
                    grid.ex(i, j, k) = grid.ex(i, j + numInternalCells.y, k);
                    grid.ey(i, j, k) = grid.ey(i, j + numInternalCells.y, k);
                    grid.ez(i, j, k) = grid.ez(i, j + numInternalCells.y, k);
                    grid.ex(i, numInternalCells.y + numGhostCells + j, k) = grid.ex(i, numGhostCells + j, k);
                    grid.ey(i, numInternalCells.y + numGhostCells + j, k) = grid.ey(i, numGhostCells + j, k);
                    grid.ez(i, numInternalCells.y + numGhostCells + j, k) = grid.ez(i, numGhostCells + j, k);
                }
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < numCells.x; i++)
            for (int j = 0; j < numCells.y; j++)
                for (int k = 0; k < numGhostCells; k++) {
                    grid.ex(i, j, k) = grid.ex(i, j, k + numInternalCells.z);
                    grid.ey(i, j, k) = grid.ey(i, j, k + numInternalCells.z);
                    grid.ez(i, j, k) = grid.ez(i, j, k + numInternalCells.z);
                    grid.ex(i, j, numInternalCells.z + numGhostCells + k) = grid.ex(i, j, numGhostCells + k);
                    grid.ey(i, j, numInternalCells.z + numGhostCells + k) = grid.ey(i, j, numGhostCells + k);
                    grid.ez(i, j, numInternalCells.z + numGhostCells + k) = grid.ez(i, j, numGhostCells + k);
                }
    }
};


#endif
