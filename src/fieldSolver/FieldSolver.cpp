#include "pica/fieldSolver/FieldSolver.h"

#include "pica/math/Constants.h"
#include "pica/grid/Grid.h"
#include "pica/threading/OpenMPHelper.h"

#include <algorithm>


namespace pica {

FieldSolver::FieldSolver(const Parameters& _parameters, Grid& _grid)
{
    parameters = _parameters;
    grid = &_grid;
}


void FieldSolver::updateA()
{
    const FP cdt = constants::c * grid->dt;
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < grid->numCells.x; i++)
    for (int j = 0; j < grid->numCells.y; j++)
    {
        #pragma simd
        for (int k = 0; k < grid->numCells.z; k++)
        {
            grid->Ax(i, j, k) -= cdt * grid->Ex(i, j, k);
            grid->Ay(i, j, k) -= cdt * grid->Ey(i, j, k);
            grid->Az(i, j, k) -= cdt * grid->Ez(i, j, k);
        }
    }
}


void FieldSolver::updateDims()
{
    updateEAreaBegin = Int3(0, 0, 0);
    updateEAreaEnd = grid->numCells - Int3(1, 1, 1);
    updateBAreaBegin = Int3(2, 2, 2);
    for (int i = parameters.dimensionality; i < 3; i++)
    {
        updateBAreaBegin[i] = 0;
        updateEAreaEnd[i]++;
    }
    updateBAreaEnd = grid->numCells - updateBAreaBegin;
}


void FieldSolver::updateInternalDims()
{
    for (int d = 0; d < 3; ++d)
    {
        internalBAreaBegin[d] = std::max(updateBAreaBegin[d], 0);
        internalBAreaEnd[d] = std::min(updateBAreaEnd[d],
            grid->numCells[d] - 0);
        internalEAreaBegin[d] = std::max(updateEAreaBegin[d], 0);
        internalEAreaEnd[d] = std::min(updateEAreaEnd[d],
            grid->numCells[d] - 0);
    }
}

} // namespace pica
