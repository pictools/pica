#include "pica/fieldSolver/Fdtd.h"

#include "pica/math/Constants.h"
#include "pica/grid/Grid.h"
#include "pica/threading/OpenMPHelper.h"

#include <algorithm>


namespace pica {

Fdtd::Fdtd(const Parameters& _parameters, Grid& _grid):
    FieldSolver(_parameters, _grid)
{
    updateDims();
    updateInternalDims();
}


// Update grid values of magnetic field in FDTD.
void Fdtd::updateHalfB()
{
    if (parameters.dimensionality == 3)
        updateHalfB3D();
    else if (parameters.dimensionality == 2)
        updateHalfB2D();
    else if (parameters.dimensionality == 1)
        updateHalfB1D();
}


void Fdtd::updateHalfB3D()
{
    updateBAreaBegin = Int3(2, 2, 2);
    updateBAreaEnd = grid->numCells - Int3(2, 2, 2);
    for (int d = 0; d < 3; ++d)
    {
        internalBAreaBegin[d] = std::max(updateBAreaBegin[d], 0);
        internalBAreaEnd[d] = std::min(updateBAreaEnd[d],
            grid->numCells[d] - 0);
    }

    const FP cdt = constants::c * grid->dt * (FP)0.5;
    const FP coeffXY = cdt / (grid->steps.x);
    const FP coeffXZ = cdt / (grid->steps.x);
    const FP coeffYX = cdt / (grid->steps.y);
    const FP coeffYZ = cdt / (grid->steps.y);
    const FP coeffZX = cdt / (grid->steps.z);
    const FP coeffZY = cdt / (grid->steps.z);

    // In central area use b(i, j, k) += c * dt * -rot(e(i, j, k)), which is:
    // b.x(i, j, k) += c * dt * ((e.y(i, j, k) - e.y(i, j, k-1)) / eps_z * dz -
    //     (e.z(i, j, k) - e.z(i, j-1, k)) / eps_y * dy),
    // b.y(i, j, k) += c * dt * ((e.z(i, j, k) - e.z(i-1, j, k)) / eps_x * dx -
    //     (e.x(i, j, k) - e.x(i, j, k-1)) / eps_z * dz),
    // b.z(i, j, k) += c * dt * ((e.x(i, j, k) - e.x(i, j-1, k)) / eps_y * dy -
    //     (e.y(i, j, k) - e.y(i-1, j, k)) / eps_x * dx),
    const Int3 begin = internalBAreaBegin;
    const Int3 end = internalBAreaEnd;
    #pragma omp parallel for collapse(2)
    for (int i = begin.x; i < end.x; i++)
    for (int j = begin.y; j < end.y; j++)
    {
        #pragma simd
        for (int k = begin.z; k < end.z; k++)
        {
            grid->Bx(i, j, k) += coeffZX * (grid->Ey(i, j, k) - grid->Ey(i, j, k - 1)) -
                coeffYX * (grid->Ez(i, j, k) - grid->Ez(i, j - 1, k));
            grid->By(i, j, k) += coeffXY * (grid->Ez(i, j, k) - grid->Ez(i - 1, j, k)) -
                coeffZY * (grid->Ex(i, j, k) - grid->Ex(i, j, k - 1));
            grid->Bz(i, j, k) += coeffYZ * (grid->Ex(i, j, k) - grid->Ex(i, j - 1, k)) -
                coeffXZ * (grid->Ey(i, j, k) - grid->Ey(i - 1, j, k));
        }
    }
}


// Update grid values of magnetic field in FDTD.
void Fdtd::updateHalfB2D()
{
    updateBAreaBegin = Int3(2, 2, 0);
    updateBAreaEnd = grid->numCells - Int3(2, 2, 0);
    for (int d = 0; d < 2; ++d)
    {
        internalBAreaBegin[d] = std::max(updateBAreaBegin[d], 0);
        internalBAreaEnd[d] = std::min(updateBAreaEnd[d],
            grid->numCells[d] - 0);
    }

    const FP cdt = constants::c * grid->dt * (FP)0.5;
    const FP coeffXY = cdt / (grid->steps.x);
    const FP coeffXZ = cdt / (grid->steps.x);
    const FP coeffYX = cdt / (grid->steps.y);
    const FP coeffYZ = cdt / (grid->steps.y);

    // In central area use b(i, j, k) += c * dt * -rot(e(i, j, k)), which is:
    // b.x(i, j, k) += c * dt * ((e.y(i, j, k) - e.y(i, j, k-1)) / eps_z * dz -
    //     (e.z(i, j, k) - e.z(i, j-1, k)) / eps_y * dy),
    // b.y(i, j, k) += c * dt * ((e.z(i, j, k) - e.z(i-1, j, k)) / eps_x * dx -
    //     (e.x(i, j, k) - e.x(i, j, k-1)) / eps_z * dz),
    // b.z(i, j, k) += c * dt * ((e.x(i, j, k) - e.x(i, j-1, k)) / eps_y * dy -
    //     (e.y(i, j, k) - e.y(i-1, j, k)) / eps_x * dx),
    const Int3 begin = internalBAreaBegin;
    const Int3 end = internalBAreaEnd;
    #pragma omp parallel for
    for (int i = begin.x; i < end.x; i++) {
        #pragma simd
        for (int j = begin.y; j < end.y; j++)
        {
            grid->Bx(i, j, 0) += -coeffYX * (grid->Ez(i, j, 0) - grid->Ez(i, j - 1, 0));
            grid->By(i, j, 0) += coeffXY * (grid->Ez(i, j, 0) - grid->Ez(i - 1, j, 0));
            grid->Bz(i, j, 0) += coeffYZ * (grid->Ex(i, j, 0) - grid->Ex(i, j - 1, 0)) -
                coeffXZ * (grid->Ey(i, j, 0) - grid->Ey(i - 1, j, 0));
        }
    }
}


// Update grid values of magnetic field in FDTD.
void Fdtd::updateHalfB1D()
{
    updateBAreaBegin = Int3(2, 0, 0);
    updateBAreaEnd = grid->numCells - Int3(2, 0, 0);
    for (int d = 0; d < 1; ++d)
    {
        internalBAreaBegin[d] = std::max(updateBAreaBegin[d], 0);
        internalBAreaEnd[d] = std::min(updateBAreaEnd[d],
            grid->numCells[d] - 0);
    }

    const FP cdt = constants::c * grid->dt * (FP)0.5;
    const FP coeffXY = cdt / (grid->steps.x);
    const FP coeffXZ = cdt / (grid->steps.x);

    // In central area use b(i, j, k) += c * dt * -rot(e(i, j, k)), which is:
    // b.x(i, j, k) += c * dt * ((e.y(i, j, k) - e.y(i, j, k-1)) / eps_z * dz -
    //     (e.z(i, j, k) - e.z(i, j-1, k)) / eps_y * dy),
    // b.y(i, j, k) += c * dt * ((e.z(i, j, k) - e.z(i-1, j, k)) / eps_x * dx -
    //     (e.x(i, j, k) - e.x(i, j, k-1)) / eps_z * dz),
    // b.z(i, j, k) += c * dt * ((e.x(i, j, k) - e.x(i, j-1, k)) / eps_y * dy -
    //     (e.y(i, j, k) - e.y(i-1, j, k)) / eps_x * dx),
    const Int3 begin = internalBAreaBegin;
    const Int3 end = internalBAreaEnd;
    #pragma omp parallel for
    for (int i = begin.x; i < end.x; i++) {
        grid->By(i, 0, 0) += coeffXY * (grid->Ez(i, 0, 0) - grid->Ez(i - 1, 0, 0));
        grid->Bz(i, 0, 0) += -coeffXZ * (grid->Ey(i, 0, 0) - grid->Ey(i - 1, 0, 0));
    }
}


// Update grid values of electric field in FDTD.
void Fdtd::updateE()
{
    if (parameters.dimensionality == 3)
        updateE3D();
    else if (parameters.dimensionality == 2)
        updateE2D();
    else if (parameters.dimensionality == 1)
        updateE1D();
}


void Fdtd::updateE3D()
{
    updateEAreaBegin = Int3(0, 0, 0);
    updateEAreaEnd = grid->numCells - Int3(1, 1, 1);
    for (int d = 0; d < 3; ++d)
    {
        internalEAreaBegin[d] = std::max(updateEAreaBegin[d], 0);
        internalEAreaEnd[d] = std::min(updateEAreaEnd[d],
            grid->numCells[d] - 0);
    }

    const FP coeffCurrent = -(FP)4 * constants::pi * grid->dt;
    const FP cdt = constants::c * grid->dt;
    const FP coeffXY = cdt / (grid->steps.x);
    const FP coeffXZ = cdt / (grid->steps.x);
    const FP coeffYX = cdt / (grid->steps.y);
    const FP coeffYZ = cdt / (grid->steps.y);
    const FP coeffZX = cdt / (grid->steps.z);
    const FP coeffZY = cdt / (grid->steps.z);

    // In internal area use:
    // e.x(i, j, k) += dt * -4pi * j.x(i, j, k) + c * dt * ((b.z(i, j+1, k) -
    //     b.z(i, j, k)) / eps_y * dy - (b.y(i, j, k+1) - b.y(i, j, k)) / eps_z * dz),
    // e.y(i, j, k) += dt * -4pi * j.y(i, j, k) + c * dt * ((b.x(i, j, k+1) -
    //     b.x(i, j, k)) / eps_z * dz - (b.z(i+1, j, k) - b.z(i, j, k)) / eps_x * dx),
    // e.z(i, j, k) += dt * -4pi * j.z(i, j, k) + c * dt * ((b.y(i+1, j, k) -
    //     b.y(i, j, k)) / eps_x * dx - (b.x(i, j+1, k) - b.x(i, j, k)) / eps_y * dy),
    const Int3 begin = internalEAreaBegin;
    const Int3 end = internalEAreaEnd;
    #pragma omp parallel for collapse(2)
    for (int i = begin.x; i < end.x; i++)
    for (int j = begin.y; j < end.y; j++)
    {
        #pragma simd
        for (int k = begin.z; k < end.z; k++)
        {
            grid->Ex(i, j, k) += coeffCurrent * grid->Jx(i, j, k) +
                coeffYX * (grid->Bz(i, j + 1, k) - grid->Bz(i, j, k)) -
                coeffZX * (grid->By(i, j, k + 1) - grid->By(i, j, k));
            grid->Ey(i, j, k) += coeffCurrent * grid->Jy(i, j, k) +
                coeffZY * (grid->Bx(i, j, k + 1) - grid->Bx(i, j, k)) -
                coeffXY * (grid->Bz(i + 1, j, k) - grid->Bz(i, j, k));
            grid->Ez(i, j, k) += coeffCurrent * grid->Jz(i, j, k) +
                coeffXZ * (grid->By(i + 1, j, k) - grid->By(i, j, k)) -
                coeffYZ * (grid->Bx(i, j + 1, k) - grid->Bx(i, j, k));
        }
    }

    // Process edge values
    if (updateEAreaEnd.x == grid->numCells.x - 1)
    {
        int i = updateEAreaEnd.x;
        #pragma omp parallel for
        for (int j = begin.y; j < end.y; j++)
        for (int k = begin.z; k < end.z; k++)
            grid->Ex(i, j, k) += coeffCurrent * grid->Jx(i, j, k) +
                coeffYX * (grid->Bz(i, j + 1, k) - grid->Bz(i, j, k)) -
                coeffZX * (grid->By(i, j, k + 1) - grid->By(i, j, k));
    }
    if (updateEAreaEnd.y == grid->numCells.y - 1)
    {
        int j = updateEAreaEnd.y;
        #pragma omp parallel for
        for (int i = begin.x; i < end.x; i++)
        for (int k = begin.z; k < end.z; k++)
            grid->Ey(i, j, k) += coeffCurrent * grid->Jy(i, j, k) +
                coeffZY * (grid->Bx(i, j, k + 1) - grid->Bx(i, j, k)) -
                coeffXY * (grid->Bz(i + 1, j, k) - grid->Bz(i, j, k));
    }
    if (updateEAreaEnd.z == grid->numCells.z - 1)
    {
        int k = updateEAreaEnd.z;
        #pragma omp parallel for
        for (int i = begin.x; i < end.x; i++)
        for (int j = begin.y; j < end.y; j++)
            grid->Ez(i, j, k) += coeffCurrent * grid->Jz(i, j, k) +
                coeffXZ * (grid->By(i + 1, j, k) - grid->By(i, j, k)) -
                coeffYZ * (grid->Bx(i, j + 1, k) - grid->Bx(i, j, k));
    }
}


void Fdtd::updateE2D()
{
    updateEAreaBegin = Int3(0, 0, 0);
    updateEAreaEnd = grid->numCells - Int3(1, 1, 0);
    for (int d = 0; d < 2; ++d)
    {
        internalEAreaBegin[d] = std::max(updateEAreaBegin[d], 0);
        internalEAreaEnd[d] = std::min(updateEAreaEnd[d],
            grid->numCells[d] - 0);
    }

    const FP coeffCurrent = -(FP)4 * constants::pi * grid->dt;
    const FP cdt = constants::c * grid->dt;
    const FP coeffXY = cdt / (grid->steps.x);
    const FP coeffXZ = cdt / (grid->steps.x);
    const FP coeffYX = cdt / (grid->steps.y);
    const FP coeffYZ = cdt / (grid->steps.y);

    // In internal area use:
    // e.x(i, j, k) += dt * -4pi * j.x(i, j, k) + c * dt * ((b.z(i, j+1, k) -
    //     b.z(i, j, k)) / eps_y * dy - (b.y(i, j, k+1) - b.y(i, j, k)) / eps_z * dz),
    // e.y(i, j, k) += dt * -4pi * j.y(i, j, k) + c * dt * ((b.x(i, j, k+1) -
    //     b.x(i, j, k)) / eps_z * dz - (b.z(i+1, j, k) - b.z(i, j, k)) / eps_x * dx),
    // e.z(i, j, k) += dt * -4pi * j.z(i, j, k) + c * dt * ((b.y(i+1, j, k) -
    //     b.y(i, j, k)) / eps_x * dx - (b.x(i, j+1, k) - b.x(i, j, k)) / eps_y * dy),
    const Int3 begin = internalEAreaBegin;
    const Int3 end = internalEAreaEnd;
    #pragma omp parallel for
    for (int i = begin.x; i < end.x; i++) {
        #pragma simd
        for (int j = begin.y; j < end.y; j++) {
            grid->Ex(i, j, 0) += coeffCurrent * grid->Jx(i, j, 0) +
                coeffYX * (grid->Bz(i, j + 1, 0) - grid->Bz(i, j, 0));
            grid->Ey(i, j, 0) += coeffCurrent * grid->Jy(i, j, 0) -
                coeffXY * (grid->Bz(i + 1, j, 0) - grid->Bz(i, j, 0));
            grid->Ez(i, j, 0) += coeffCurrent * grid->Jz(i, j, 0) +
                coeffXZ * (grid->By(i + 1, j, 0) - grid->By(i, j, 0)) -
                coeffYZ * (grid->Bx(i, j + 1, 0) - grid->Bx(i, j, 0));
        }
    }

    // Process edge values
    if (updateEAreaEnd.x == grid->numCells.x - 1)
    {
        int i = updateEAreaEnd.x;
        #pragma omp parallel for
        for (int j = begin.y; j < end.y; j++)
            grid->Ex(i, j, 0) += coeffCurrent * grid->Jx(i, j, 0) +
                coeffYX * (grid->Bz(i, j + 1, 0) - grid->Bz(i, j, 0));
    }
    if (updateEAreaEnd.y == grid->numCells.y - 1)
    {
        int j = updateEAreaEnd.y;
        #pragma omp parallel for
        for (int i = begin.x; i < end.x; i++)
            grid->Ey(i, j, 0) += coeffCurrent * grid->Jy(i, j, 0) -
                coeffXY * (grid->Bz(i + 1, j, 0) - grid->Bz(i, j, 0));
    }
}


void Fdtd::updateE1D()
{
    updateEAreaBegin = Int3(0, 0, 0);
    updateEAreaEnd = grid->numCells - Int3(1, 0, 0);
    for (int d = 0; d < 1; ++d)
    {
        internalEAreaBegin[d] = std::max(updateEAreaBegin[d], 0);
        internalEAreaEnd[d] = std::min(updateEAreaEnd[d],
            grid->numCells[d] - 0);
    }

    const FP coeffCurrent = -(FP)4 * constants::pi * grid->dt;
    const FP cdt = constants::c * grid->dt;
    const FP coeffXY = cdt / (grid->steps.x);
    const FP coeffXZ = cdt / (grid->steps.x);

    // In internal area use:
    // e.x(i, j, k) += dt * -4pi * j.x(i, j, k) + c * dt * ((b.z(i, j+1, k) -
    //     b.z(i, j, k)) / eps_y * dy - (b.y(i, j, k+1) - b.y(i, j, k)) / eps_z * dz),
    // e.y(i, j, k) += dt * -4pi * j.y(i, j, k) + c * dt * ((b.x(i, j, k+1) -
    //     b.x(i, j, k)) / eps_z * dz - (b.z(i+1, j, k) - b.z(i, j, k)) / eps_x * dx),
    // e.z(i, j, k) += dt * -4pi * j.z(i, j, k) + c * dt * ((b.y(i+1, j, k) -
    //     b.y(i, j, k)) / eps_x * dx - (b.x(i, j+1, k) - b.x(i, j, k)) / eps_y * dy),
    const Int3 begin = internalEAreaBegin;
    const Int3 end = internalEAreaEnd;
    #pragma omp parallel for
    for (int i = begin.x; i < end.x; i++) {
        grid->Ex(i, 0, 0) += coeffCurrent * grid->Jx(i, 0, 0);
        grid->Ey(i, 0, 0) += coeffCurrent * grid->Jy(i, 0, 0) -
            coeffXY * (grid->Bz(i + 1, 0, 0) - grid->Bz(i, 0, 0));
        grid->Ez(i, 0, 0) += coeffCurrent * grid->Jz(i, 0, 0) +
            coeffXZ * (grid->By(i + 1, 0, 0) - grid->By(i, 0, 0));
    }
}

} // namespace pica
