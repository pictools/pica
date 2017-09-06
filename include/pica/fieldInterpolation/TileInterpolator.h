#ifndef PICA_TILEINTERPOLATOR_H
#define PICA_TILEINTERPOLATOR_H


#include "pica/fieldInterpolation/Formfactor.h"
#include "pica/grid/Grid.h"
#include "pica/threading/OpenMPHelper.h"


namespace pica {

/* Base class for tile field interpolators. It is responsible for allocating
memory and pre-loading field values used for particles of a tile. */
template<int TileSize>
class TileInterpolator
{
public:

    static const int tileSize = TileSize;

    TileInterpolator(int tileI, int tileJ, int tileK, const Grid* _grid)
        : grid(_grid)
    {
        int tileOffset = (tileSize + 1) / 2;
        Int3 baseGridIdx = Int3(tileI, tileJ, tileK) -
            Int3(tileOffset, tileOffset, tileOffset)
            + grid->getNumExternalLeftCells();
        for (int i = 0; i < tileSize; ++i)
        for (int j = 0; j < tileSize; ++j)
        for (int k = 0; k < tileSize; ++k)
        {
            Int3 gridIdx = remainder(Int3(i, j, k) + baseGridIdx, grid->numCells);
            ex[i][j][k] = grid->Ex(gridIdx);
            ey[i][j][k] = grid->Ey(gridIdx);
            ez[i][j][k] = grid->Ez(gridIdx);
            bx[i][j][k] = grid->Bx(gridIdx);
            by[i][j][k] = grid->By(gridIdx);
            bz[i][j][k] = grid->Bz(gridIdx);
        }
        invSteps = FP3(1.0, 1.0, 1.0) / grid->steps;
        normalizedCellOrigin = grid->origin / grid->steps + FP3(baseGridIdx);
        normalizedMiddleCellOrigin = normalizedCellOrigin + FP3(0.5, 0.5, 0.5);
    }

    // Default implementation is to just call the corresponding of the grid
    void getFields(const FP3& coords, FP3& e, FP3& b) const
    {
        grid->getFields(coords.x, coords.y, coords.z, e, b);
    }

protected:

    const Grid *grid;
    FP ex[tileSize][tileSize][tileSize];
    FP ey[tileSize][tileSize][tileSize];
    FP ez[tileSize][tileSize][tileSize];
    FP bx[tileSize][tileSize][tileSize];
    FP by[tileSize][tileSize][tileSize];
    FP bz[tileSize][tileSize][tileSize];
    FP3 invSteps; // 1 / steps
    FP3 normalizedCellOrigin, normalizedMiddleCellOrigin;

};


class TileInterpolatorCIC: public TileInterpolator<3>
{
public:

    TileInterpolatorCIC(int tileI, int tileJ, int tileK, const Grid* _grid):
        TileInterpolator<3>(tileI, tileJ, tileK, _grid) {}

    void getFields(const FP3& coords, FP3& e, FP3& b) const
    {
        FP3 normalizedCoords = coords * invSteps;
        FP3 internalCoords = normalizedCoords - normalizedCellOrigin;
        Int3 gridIdx = truncate(internalCoords);
        FP3 coeff = internalCoords - FP3(gridIdx);
        internalCoords = normalizedCoords - normalizedMiddleCellOrigin;
        Int3 middleGridIdx = truncate(internalCoords);
        FP3 middleCoeff = internalCoords - FP3(middleGridIdx);
        e.x = getFieldComponent(coeff.x, middleCoeff.y, middleCoeff.z,
            gridIdx.x, middleGridIdx.y, middleGridIdx.z, ex);
        e.y = getFieldComponent(middleCoeff.x, coeff.y, middleCoeff.z,
            middleGridIdx.x, gridIdx.y, middleGridIdx.z, ey);
        e.z = getFieldComponent(middleCoeff.x, middleCoeff.y, coeff.z,
            middleGridIdx.x, middleGridIdx.y, gridIdx.z, ez);
        b.x = getFieldComponent(middleCoeff.x, coeff.y, coeff.z,
            middleGridIdx.x, gridIdx.y, gridIdx.z, bx);
        b.y = getFieldComponent(coeff.x, middleCoeff.y, coeff.z,
            gridIdx.x, middleGridIdx.y, gridIdx.z, by);
        b.z = getFieldComponent(coeff.x, coeff.y, middleCoeff.z,
            gridIdx.x, gridIdx.y, middleGridIdx.z, bz);
    }

private:

    // Interpolate field f from base index (i, j, k) and coefficients (x, y, z).
    FP getFieldComponent(FP x, FP y, FP z, int i, int j, int k,
        const FP f[tileSize][tileSize][tileSize]) const
    {
        const FP ix = (FP)1 - x;
        const FP iy = (FP)1 - y;
        const FP iz = (FP)1 - z;
        return
            ix * (iy * (iz * f[i][j][k] + z * f[i][j][k + 1]) +
                   y * (iz * f[i][j + 1][k] + z * f[i][j + 1][k + 1])) +
             x * (iy * (iz * f[i + 1][j][k] + z * f[i + 1][j][k + 1]) +
                   y * (iz * f[i + 1][j + 1][k] + z * f[i + 1][j + 1][k + 1]));
    }

};


class TileInterpolatorTSC: public TileInterpolator<4>
{
public:

    TileInterpolatorTSC(int tileI, int tileJ, int tileK, const Grid* _grid):
        TileInterpolator<4>(tileI, tileJ, tileK, _grid) {}

    void getFields(const FP3& coords, FP3& e, FP3& b) const
    {
        FP3 normalizedCoords = coords * invSteps;
        FP3 internalCoords = normalizedCoords - normalizedCellOrigin;
        Int3 gridIdx = truncate(internalCoords + FP3(0.5, 0.5, 0.5));
        FP3 coeff = internalCoords - FP3(gridIdx);
        internalCoords = normalizedCoords - normalizedMiddleCellOrigin;
        Int3 middleGridIdx = truncate(internalCoords + FP3(0.5, 0.5, 0.5));
        FP3 middleCoeff = internalCoords - FP3(middleGridIdx);

        e.x = getFieldComponent(coeff.x, middleCoeff.y, middleCoeff.z,
            gridIdx.x, middleGridIdx.y, middleGridIdx.z, ex);
        e.y = getFieldComponent(middleCoeff.x, coeff.y, middleCoeff.z,
            middleGridIdx.x, gridIdx.y, middleGridIdx.z, ey);
        e.z = getFieldComponent(middleCoeff.x, middleCoeff.y, coeff.z,
            middleGridIdx.x, middleGridIdx.y, gridIdx.z, ez);
        b.x = getFieldComponent(middleCoeff.x, coeff.y, coeff.z,
            middleGridIdx.x, gridIdx.y, gridIdx.z, bx);
        b.y = getFieldComponent(coeff.x, middleCoeff.y, coeff.z,
            gridIdx.x, middleGridIdx.y, gridIdx.z, by);
        b.z = getFieldComponent(coeff.x, coeff.y, middleCoeff.z,
            gridIdx.x, gridIdx.y, middleGridIdx.z, bz);
    }

private:

    FP getFieldComponent(FP x, FP y, FP z, int i, int j, int k,
        const FP f[tileSize][tileSize][tileSize]) const
    {
        FP cx[3], cy[3], cz[3];
        for (int ii = 0; ii < 3; ii++)
            cx[ii] = formfactorTSC(FP(ii - 1) - x);
        for (int jj = 0; jj < 3; jj++)
            cy[jj] = formfactorTSC(FP(jj - 1) - y);
        for (int kk = 0; kk < 3; kk++)
            cz[kk] = formfactorTSC(FP(kk - 1) - z);
        FP result = 0;
        for (int ii = -1; ii <= 1; ii++)
        for (int jj = -1; jj <= 1; jj++)
        for (int kk = -1; kk <= 1; kk++)
            result += cx[ii + 1] * cy[jj + 1] * cz[kk + 1] * f[i + ii][j + jj][k + kk];
        return result;
    }

};


/* Base class for tile field interpolators in XY 2D case. It is responsible for allocating
memory and pre-loading field values used for particles of a tile. */
template<int TileSize>
class XYTileInterpolator
{
public:
    
    static const int tileSize = TileSize;

    XYTileInterpolator(int tileI, int tileJ, int tileK, const Grid* _grid)
        : grid(_grid)
    {
        int tileOffset = (tileSize + 1) / 2;
        Int3 baseGridIdx = Int3(tileI, tileJ, 0) -
            Int3(tileOffset, tileOffset, 0)
            + grid->getNumExternalLeftCells();
        for (int i = 0; i < tileSize; ++i)
        for (int j = 0; j < tileSize; ++j)
        {
            Int3 gridIdx = Int3(i, j, 0) + baseGridIdx;
            ex[i][j] = grid->Ex(gridIdx);
            ey[i][j] = grid->Ey(gridIdx);
            ez[i][j] = grid->Ez(gridIdx);
            bx[i][j] = grid->Bx(gridIdx);
            by[i][j] = grid->By(gridIdx);
            bz[i][j] = grid->Bz(gridIdx);
        }
        invSteps = FP3(1.0, 1.0, 1.0) / grid->steps;
        normalizedCellOrigin = grid->origin / grid->steps + FP3(baseGridIdx);
        normalizedMiddleCellOrigin = normalizedCellOrigin + FP3(0.5, 0.5, 0.5);
    }

    // Default implementation is to just call the corresponding of the grid
    void getFields(const FP3& coords, FP3& e, FP3& b) const
    {
        grid->getFields(coords.x, coords.y, coords.z, e, b);
    }

protected:

    const Grid *grid;
    FP ex[tileSize][tileSize];
    FP ey[tileSize][tileSize];
    FP ez[tileSize][tileSize];
    FP bx[tileSize][tileSize];
    FP by[tileSize][tileSize];
    FP bz[tileSize][tileSize];
    FP3 invSteps; // 1 / steps
    FP3 normalizedCellOrigin, normalizedMiddleCellOrigin;

};


class XYTileInterpolatorTSC: public XYTileInterpolator<4>
{
public:

    XYTileInterpolatorTSC(int tileI, int tileJ, int tileK, const Grid* _grid):
        XYTileInterpolator<4>(tileI, tileJ, tileK, _grid) {}

    void getFields(const FP3& coords, FP3& e, FP3& b) const
    {
        FP3 normalizedCoords(coords.x * invSteps.x, coords.y * invSteps.y, 0);
        FP3 internalCoords(normalizedCoords.x - normalizedCellOrigin.x,
                           normalizedCoords.y - normalizedCellOrigin.y, 0);
        Int3 gridIdx((int)(internalCoords.x + (FP)0.5),
                     (int)(internalCoords.y + (FP)0.5), 0);
        FP3 coeff(internalCoords.x - (FP)gridIdx.x,
                  internalCoords.y - (FP)gridIdx.y, 0);
        internalCoords.x = normalizedCoords.x - normalizedMiddleCellOrigin.x;
        internalCoords.y = normalizedCoords.y - normalizedMiddleCellOrigin.y;
        Int3 middleGridIdx((int)(internalCoords.x + (FP)0.5),
                           (int)(internalCoords.y + (FP)0.5), 0);
        FP3 middleCoeff(internalCoords.x - (FP)middleGridIdx.x,
                        internalCoords.y - (FP)middleGridIdx.y, 0);

        e.x = getFieldComponent(coeff.x, middleCoeff.y, 
            gridIdx.x, middleGridIdx.y, ex);
        e.y = getFieldComponent(middleCoeff.x, coeff.y,
            middleGridIdx.x, gridIdx.y, ey);
        e.z = getFieldComponent(middleCoeff.x, middleCoeff.y,
            middleGridIdx.x, middleGridIdx.y, ez);
        b.x = getFieldComponent(middleCoeff.x, coeff.y,
            middleGridIdx.x, gridIdx.y, bx);
        b.y = getFieldComponent(coeff.x, middleCoeff.y,
            gridIdx.x, middleGridIdx.y, by);
        b.z = getFieldComponent(coeff.x, coeff.y,
            gridIdx.x, gridIdx.y, bz);
    }

private:

    FP getFieldComponent(FP x, FP y, int i, int j,
        const FP f[tileSize][tileSize]) const
    {
        FP cx[3], cy[3];
        for (int ii = 0; ii < 3; ii++)
            cx[ii] = formfactorTSC(FP(ii - 1) - x);
        for (int jj = 0; jj < 3; jj++)
            cy[jj] = formfactorTSC(FP(jj - 1) - y);
        FP result = 0;
        for (int ii = -1; ii <= 1; ii++)
        for (int jj = -1; jj <= 1; jj++)
            result += cx[ii + 1] * cy[jj + 1] * f[i + ii][j + jj];
        return result;
    }

};


/* Base class for tile field interpolators in X 1D case. It is responsible for allocating
memory and pre-loading field values used for particles of a tile. */
template<int TileSize>
class XTileInterpolator
{
public:
    
    static const int tileSize = TileSize;

    XTileInterpolator(int tileI, int tileJ, int tileK, const Grid* _grid)
        : grid(_grid)
    {
        int tileOffset = (tileSize + 1) / 2;
        Int3 baseGridIdx = Int3(tileI, 0, 0) - Int3(tileOffset, 0, 0)
            + grid->getNumExternalLeftCells();
        for (int i = 0; i < tileSize; ++i)
        {
            Int3 gridIdx = Int3(i, 0, 0) + baseGridIdx;
            ex[i] = grid->Ex(gridIdx);
            ey[i] = grid->Ey(gridIdx);
            ez[i] = grid->Ez(gridIdx);
            bx[i] = grid->Bx(gridIdx);
            by[i] = grid->By(gridIdx);
            bz[i] = grid->Bz(gridIdx);
        }
        invStep = (FP)1.0 / grid->steps.x;
        normalizedCellOrigin = grid->origin.x / grid->steps.x + FP(baseGridIdx.x);
        normalizedMiddleCellOrigin = normalizedCellOrigin + (FP)0.5;
    }

    // Default implementation is to just call the corresponding of the grid
    void getFields(const FP3& coords, FP3& e, FP3& b) const
    {
        grid->getFields(coords.x, coords.y, coords.z, e, b);
    }

protected:

    const Grid *grid;
    FP ex[tileSize];
    FP ey[tileSize];
    FP ez[tileSize];
    FP bx[tileSize];
    FP by[tileSize];
    FP bz[tileSize];
    FP invStep; // 1 / steps.x
    FP normalizedCellOrigin, normalizedMiddleCellOrigin;

};


class XTileInterpolatorTSC: public XTileInterpolator<4>
{
public:

    XTileInterpolatorTSC(int tileI, int tileJ, int tileK, const Grid* _grid):
        XTileInterpolator<4>(tileI, tileJ, tileK, _grid) {}

    void getFields(const FP3& coords, FP3& e, FP3& b) const
    {
        FP normalizedCoords = coords.x * invStep;
        FP internalCoords = normalizedCoords - normalizedCellOrigin;
        int gridIdx = (int)(internalCoords + (FP)0.5);
        FP coeff = internalCoords - (FP)gridIdx;
        internalCoords = normalizedCoords - normalizedMiddleCellOrigin;
        int middleGridIdx = (int)(internalCoords + (FP)0.5);
        FP middleCoeff = internalCoords - (FP)middleGridIdx;

        e.x = getFieldComponent(coeff, gridIdx, ex);
        e.y = getFieldComponent(middleCoeff, middleGridIdx, ey);
        e.z = getFieldComponent(middleCoeff, middleGridIdx, ez);
        b.x = getFieldComponent(middleCoeff, middleGridIdx, bx);
        b.y = getFieldComponent(coeff, gridIdx, by);
        b.z = getFieldComponent(coeff, gridIdx, bz);
    }

private:

    FP getFieldComponent(FP x, int i, const FP f[tileSize]) const
    {
        return formfactorTSC((FP)1.0 + x) * f[i - 1] +
            formfactorTSC(x) * f[i] + formfactorTSC((FP)1.0 - x) * f[i + 1];
    }

};

} // namespace pica


#endif
