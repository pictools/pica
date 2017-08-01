#include "pica/Grid.h"
#include "pica/OpenMPHelper.h"


namespace pica {

Grid::Grid(const Int3 & _numInternalCells, FP _dt, const FP3 & minCoords, const FP3 & _steps, 
           const Int3 & _globalGridDims, const Int3 & _domainOffset):
    globalGridDims(_globalGridDims),
    domainOffset(_domainOffset),
    steps(_steps),
    dt(_dt),
    numInternalCells(_numInternalCells),
    numCells(numInternalCells + getNumExternalLeftCells() + getNumExternalRightCells()),
    Ex(numCells), Ey(numCells), Ez(numCells),
    Bx(numCells), By(numCells), Bz(numCells),
    Jx(numCells), Jy(numCells), Jz(numCells),
    Ax(numCells), Ay(numCells), Az(numCells),
    shiftEJx(FP3(0, 0.5, 0.5) * steps),
    shiftEJy(FP3(0.5, 0, 0.5) * steps),
    shiftEJz(FP3(0.5, 0.5, 0) * steps),
    shiftBx(FP3(0.5, 0, 0) * steps),
    shiftBy(FP3(0, 0.5, 0) * steps),
    shiftBz(FP3(0, 0, 0.5) * steps),
    shiftRho(steps * 0.5),
    origin(minCoords - steps * getNumExternalLeftCells()),
    EB_timeShift(dt / 2)
{
    setInterpolationType(Interpolation_CIC);
    setDepositionType(Deposition_CIC);
}


void Grid::getFields(FP x, FP y, FP z, FP3 & e, FP3 & b) const
{
    (this->*interpolationFunction)(x, y, z, e, b);
}


void Grid::getFieldsCIC(FP x, FP y, FP z, FP3 & e, FP3 & b) const
{
    /* For each component of E and B get grid index and internal coords,
       use it as base index and coefficients of interpolation. */
    FP3 coords(x, y, z);
    Int3 idx;
    FP3 internalCoords;

    getGridCoords(coords, shiftEJx, idx, internalCoords);
    e.x = Ex.interpolateCIC(idx, internalCoords);
    getGridCoords(coords, shiftEJy, idx, internalCoords);
    e.y = Ey.interpolateCIC(idx, internalCoords);
    getGridCoords(coords, shiftEJz, idx, internalCoords);
    e.z = Ez.interpolateCIC(idx, internalCoords);

    getGridCoords(coords, shiftBx, idx, internalCoords); 
    b.x = Bx.interpolateCIC(idx, internalCoords);
    getGridCoords(coords, shiftBy, idx, internalCoords);
    b.y = By.interpolateCIC(idx, internalCoords);
    getGridCoords(coords, shiftBz, idx, internalCoords);
    b.z = Bz.interpolateCIC(idx, internalCoords);
}


void Grid::getFieldsTSC(FP x, FP y, FP z, FP3 & e, FP3 & b) const
{
    FP3 coords(x, y, z);
    Int3 idx;
    FP3 internalCoords;

    getClosestGridCoords(coords, shiftEJx, idx, internalCoords);
    e.x = Ex.interpolateTSC(idx, internalCoords);
    getClosestGridCoords(coords, shiftEJy, idx, internalCoords);
    e.y = Ey.interpolateTSC(idx, internalCoords);
    getClosestGridCoords(coords, shiftEJz, idx, internalCoords);
    e.z = Ez.interpolateTSC(idx, internalCoords);

    getClosestGridCoords(coords, shiftBx, idx, internalCoords); 
    b.x = Bx.interpolateTSC(idx, internalCoords);
    getClosestGridCoords(coords, shiftBy, idx, internalCoords);
    b.y = By.interpolateTSC(idx, internalCoords);
    getClosestGridCoords(coords, shiftBz, idx, internalCoords);
    b.z = Bz.interpolateTSC(idx, internalCoords);
}


void Grid::getFieldsSecondOrder(FP x, FP y, FP z, FP3 & e, FP3 & b) const
{
    FP3 coords(x, y, z);
    Int3 idx;
    FP3 internalCoords;

    getClosestGridCoords(coords, shiftEJx, idx, internalCoords);
    e.x = Ex.interpolateSecondOrder(idx, internalCoords);
    getClosestGridCoords(coords, shiftEJy, idx, internalCoords);
    e.y = Ey.interpolateSecondOrder(idx, internalCoords);
    getClosestGridCoords(coords, shiftEJz, idx, internalCoords);
    e.z = Ez.interpolateSecondOrder(idx, internalCoords);

    getClosestGridCoords(coords, shiftBx, idx, internalCoords); 
    b.x = Bx.interpolateSecondOrder(idx, internalCoords);
    getClosestGridCoords(coords, shiftBy, idx, internalCoords);
    b.y = By.interpolateSecondOrder(idx, internalCoords);
    getClosestGridCoords(coords, shiftBz, idx, internalCoords);
    b.z = Bz.interpolateSecondOrder(idx, internalCoords);
}


void Grid::getFieldsFourthOrder(FP x, FP y, FP z, FP3 & e, FP3 & b) const
{
    FP3 coords(x, y, z);
    Int3 idx;
    FP3 internalCoords;

    getClosestGridCoords(coords, shiftEJx, idx, internalCoords);
    e.x = Ex.interpolateFourthOrder(idx, internalCoords);
    getClosestGridCoords(coords, shiftEJy, idx, internalCoords);
    e.y = Ey.interpolateFourthOrder(idx, internalCoords);
    getClosestGridCoords(coords, shiftEJz, idx, internalCoords);
    e.z = Ez.interpolateFourthOrder(idx, internalCoords);

    getClosestGridCoords(coords, shiftBx, idx, internalCoords); 
    b.x = Bx.interpolateFourthOrder(idx, internalCoords);
    getClosestGridCoords(coords, shiftBy, idx, internalCoords);
    b.y = By.interpolateFourthOrder(idx, internalCoords);
    getClosestGridCoords(coords, shiftBz, idx, internalCoords);
    b.z = Bz.interpolateFourthOrder(idx, internalCoords);
}


void Grid::getFieldsPCS(FP x, FP y, FP z, FP3 & e, FP3 & b) const
{
    FP3 coords(x, y, z);
    Int3 idx;
    FP3 internalCoords;

    getGridCoords(coords, shiftEJx, idx, internalCoords);
    e.x = Ex.interpolatePCS(idx, internalCoords);
    getGridCoords(coords, shiftEJy, idx, internalCoords);
    e.y = Ey.interpolatePCS(idx, internalCoords);
    getGridCoords(coords, shiftEJz, idx, internalCoords);
    e.z = Ez.interpolatePCS(idx, internalCoords);

    getGridCoords(coords, shiftBx, idx, internalCoords); 
    b.x = Bx.interpolatePCS(idx, internalCoords);
    getGridCoords(coords, shiftBy, idx, internalCoords);
    b.y = By.interpolatePCS(idx, internalCoords);
    getGridCoords(coords, shiftBz, idx, internalCoords);
    b.z = Bz.interpolatePCS(idx, internalCoords);
}


void Grid::getJ(FP x, FP y, FP z, FP3 & j) const
{
    /* For each component of J get grid index and internal coords,
       use it as base index and coefficients of interpolation. */
    FP3 coords(x, y, z);
    Int3 idx;
    FP3 internalCoords;

    getGridCoords(coords, shiftEJx, idx, internalCoords);
    j.x = Jx.interpolateCIC(idx, internalCoords);
    getGridCoords(coords, shiftEJy, idx, internalCoords);
    j.y = Jy.interpolateCIC(idx, internalCoords);
    getGridCoords(coords, shiftEJz, idx, internalCoords);
    j.z = Jz.interpolateCIC(idx, internalCoords);
}


void Grid::getA(FP x, FP y, FP z, FP3 & a) const
{
    FP3 coords(x, y, z);
    Int3 idx;
    FP3 internalCoords;

    getGridCoords(coords, shiftEJx, idx, internalCoords);
    a.x = Ax.interpolateCIC(idx, internalCoords);
    getGridCoords(coords, shiftEJy, idx, internalCoords);
    a.y = Ay.interpolateCIC(idx, internalCoords);
    getGridCoords(coords, shiftEJz, idx, internalCoords);
    a.z = Az.interpolateCIC(idx, internalCoords);
}


FP Grid::getExCIC(const FP3& coords) const
{
    Int3 idx;
    FP3 internalCoords;
    getGridCoords(coords, shiftEJx, idx, internalCoords);
    return Ex.interpolateCIC(idx, internalCoords);
}


FP Grid::getEyCIC(const FP3& coords) const
{
    Int3 idx;
    FP3 internalCoords;
    getGridCoords(coords, shiftEJy, idx, internalCoords);
    return Ey.interpolateCIC(idx, internalCoords);
}


FP Grid::getEzCIC(const FP3& coords) const
{
    Int3 idx;
    FP3 internalCoords;
    getGridCoords(coords, shiftEJz, idx, internalCoords);
    return Ez.interpolateCIC(idx, internalCoords);
}


FP Grid::getBxCIC(const FP3& coords) const
{
    Int3 idx;
    FP3 internalCoords;
    getGridCoords(coords, shiftBx, idx, internalCoords);
    return Bx.interpolateCIC(idx, internalCoords);
}


FP Grid::getByCIC(const FP3& coords) const
{
    Int3 idx;
    FP3 internalCoords;
    getGridCoords(coords, shiftBy, idx, internalCoords);
    return By.interpolateCIC(idx, internalCoords);
}


FP Grid::getBzCIC(const FP3& coords) const
{
    Int3 idx;
    FP3 internalCoords;
    getGridCoords(coords, shiftBz, idx, internalCoords);
    return Bz.interpolateCIC(idx, internalCoords);
}


FP Grid::getJxCIC(const FP3& coords) const
{
    Int3 idx;
    FP3 internalCoords;
    getGridCoords(coords, shiftEJx, idx, internalCoords);
    return Jx.interpolateCIC(idx, internalCoords);
}


FP Grid::getJyCIC(const FP3& coords) const
{
    Int3 idx;
    FP3 internalCoords;
    getGridCoords(coords, shiftEJy, idx, internalCoords);
    return Jy.interpolateCIC(idx, internalCoords);
}


FP Grid::getJzCIC(const FP3& coords) const
{
    Int3 idx;
    FP3 internalCoords;
    getGridCoords(coords, shiftEJz, idx, internalCoords);
    return Jz.interpolateCIC(idx, internalCoords);
}


FP Grid::getAxCIC(const FP3& coords) const
{
    Int3 idx;
    FP3 internalCoords;
    getGridCoords(coords, shiftEJx, idx, internalCoords);
    return Ax.interpolateCIC(idx, internalCoords);
}


FP Grid::getAyCIC(const FP3& coords) const
{
    Int3 idx;
    FP3 internalCoords;
    getGridCoords(coords, shiftEJy, idx, internalCoords);
    return Ay.interpolateCIC(idx, internalCoords);
}


FP Grid::getAzCIC(const FP3& coords) const
{
    Int3 idx;
    FP3 internalCoords;
    getGridCoords(coords, shiftEJz, idx, internalCoords);
    return Az.interpolateCIC(idx, internalCoords);
}


FP Grid::getExSecondOrder(const FP3& coords) const
{
    Int3 idx;
    FP3 internalCoords;
    getClosestGridCoords(coords, shiftEJx, idx, internalCoords);
    return Ex.interpolateSecondOrder(idx, internalCoords);
}


FP Grid::getEySecondOrder(const FP3& coords) const
{
    Int3 idx;
    FP3 internalCoords;
    getClosestGridCoords(coords, shiftEJy, idx, internalCoords);
    return Ey.interpolateSecondOrder(idx, internalCoords);
}


FP Grid::getEzSecondOrder(const FP3& coords) const
{
    Int3 idx;
    FP3 internalCoords;
    getClosestGridCoords(coords, shiftEJz, idx, internalCoords);
    return Ez.interpolateSecondOrder(idx, internalCoords);
}


FP Grid::getBxSecondOrder(const FP3& coords) const
{
    Int3 idx;
    FP3 internalCoords;
    getClosestGridCoords(coords, shiftBx, idx, internalCoords);
    return Bx.interpolateSecondOrder(idx, internalCoords);
}


FP Grid::getBySecondOrder(const FP3& coords) const
{
    Int3 idx;
    FP3 internalCoords;
    getClosestGridCoords(coords, shiftBy, idx, internalCoords);
    return By.interpolateSecondOrder(idx, internalCoords);
}


FP Grid::getBzSecondOrder(const FP3& coords) const
{
    Int3 idx;
    FP3 internalCoords;
    getClosestGridCoords(coords, shiftBz, idx, internalCoords);
    return Bz.interpolateSecondOrder(idx, internalCoords);
}


FP Grid::getJxSecondOrder(const FP3& coords) const
{
    Int3 idx;
    FP3 internalCoords;
    getClosestGridCoords(coords, shiftEJx, idx, internalCoords);
    return Jx.interpolateSecondOrder(idx, internalCoords);
}


FP Grid::getJySecondOrder(const FP3& coords) const
{
    Int3 idx;
    FP3 internalCoords;
    getClosestGridCoords(coords, shiftEJy, idx, internalCoords);
    return Jy.interpolateSecondOrder(idx, internalCoords);
}


FP Grid::getJzSecondOrder(const FP3& coords) const
{
    Int3 idx;
    FP3 internalCoords;
    getClosestGridCoords(coords, shiftEJz, idx, internalCoords);
    return Jz.interpolateSecondOrder(idx, internalCoords);
}


FP Grid::getAxSecondOrder(const FP3& coords) const
{
    Int3 idx;
    FP3 internalCoords;
    getClosestGridCoords(coords, shiftEJx, idx, internalCoords);
    return Ax.interpolateSecondOrder(idx, internalCoords);
}


FP Grid::getAySecondOrder(const FP3& coords) const
{
    Int3 idx;
    FP3 internalCoords;
    getClosestGridCoords(coords, shiftEJy, idx, internalCoords);
    return Ay.interpolateSecondOrder(idx, internalCoords);
}


FP Grid::getAzSecondOrder(const FP3& coords) const
{
    Int3 idx;
    FP3 internalCoords;
    getClosestGridCoords(coords, shiftEJz, idx, internalCoords);
    return Az.interpolateSecondOrder(idx, internalCoords);
}


void Grid::depositJ(const FP3 & value, const FP3 & coords)
{
    (this->*depositionFunction)(value, coords);
}


void Grid::depositJCIC(const FP3 & value, const FP3 & coords)
{
    Int3 idx;
    FP3 internalCoords;
    getGridCoords(coords, shiftEJx, idx, internalCoords);
    Jx.depositCIC(value.x, idx, internalCoords);
    getGridCoords(coords, shiftEJy, idx, internalCoords);
    Jy.depositCIC(value.y, idx, internalCoords);
    getGridCoords(coords, shiftEJz, idx, internalCoords);
    Jz.depositCIC(value.z, idx, internalCoords);
}


void Grid::depositJTSC(const FP3 & value, const FP3 & coords)
{
    Int3 idx;
    FP3 internalCoords;
    getClosestGridCoords(coords, shiftEJx, idx, internalCoords);
    Jx.depositTSC(value.x, idx, internalCoords);
    getClosestGridCoords(coords, shiftEJy, idx, internalCoords);
    Jy.depositTSC(value.y, idx, internalCoords);
    getClosestGridCoords(coords, shiftEJz, idx, internalCoords);
    Jz.depositTSC(value.z, idx, internalCoords);
}


void Grid::depositJPCS(const FP3 & value, const FP3 & coords)
{
    Int3 idx;
    FP3 internalCoords;
    getGridCoords(coords, shiftEJx, idx, internalCoords);
    Jx.depositPCS(value.x, idx, internalCoords);
    getGridCoords(coords, shiftEJy, idx, internalCoords);
    Jy.depositPCS(value.y, idx, internalCoords);
    getGridCoords(coords, shiftEJz, idx, internalCoords);
    Jz.depositPCS(value.z, idx, internalCoords);
}


void Grid::zeroizeJ()
{
    Jx.zeroize();
    Jy.zeroize();
    Jz.zeroize();
}


void Grid::setInterpolationType(Grid::InterpolationType type)
{
    interpolationType = type;
    switch (interpolationType)
    {
        case Grid::Interpolation_CIC:
            interpolationFunction = &Grid::getFieldsCIC; break;
        case Grid::Interpolation_TSC:
            interpolationFunction = &Grid::getFieldsTSC; break;
        case Grid::Interpolation_PCS:
            interpolationFunction = &Grid::getFieldsPCS; break;
        case Grid::Interpolation_SecondOrder:
            interpolationFunction = &Grid::getFieldsSecondOrder; break;
        case Grid::Interpolation_FourthOrder:
            interpolationFunction = &Grid::getFieldsFourthOrder; break;
    }
}


Grid::InterpolationType Grid::getInterpolationType() const
{
    return interpolationType;
}


void Grid::setDepositionType(Grid::DepositionType type)
{
    depositionType = type;
    switch (depositionType)
    {
        case Grid::Deposition_CIC:
        case Grid::Deposition_VillasenorBuneman:
        case Grid::Deposition_ZigzagFirstOrder:
            depositionFunction = &Grid::depositJCIC; break;
        case Grid::Deposition_TSC:
        case Grid::Deposition_ZigzagSecondOrder:
        case Grid::Deposition_Esirkepov:
            depositionFunction = &Grid::depositJTSC; break;
    }
}


Grid::DepositionType Grid::getDepositionType() const
{
    return depositionType;
}


void Grid::addCurrents(const FP3 * currents, const Int3 * minCellIdx, const Int3 * maxCellIdx)
{
    Int3 numCells = *maxCellIdx - *minCellIdx;
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < numCells.x; ++i)
    for (int j = 0; j < numCells.y; ++j)
    for (int k = 0; k < numCells.z; ++k)
    {
        int idx = numCells.y * numCells.z * i + numCells.z * j + k;
        Int3 nodeIdx = *minCellIdx + Int3(i, j, k);
        Jx(nodeIdx) += currents[idx].x;
        Jy(nodeIdx) += currents[idx].y;
        Jz(nodeIdx) += currents[idx].z;
    }
}

void Grid::dumpB(FP3 * b, const Int3 * minCellIdx, const Int3 * maxCellIdx)
    const
{
    Int3 numCells = *maxCellIdx - *minCellIdx;
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < numCells.x; ++i)
    for (int j = 0; j < numCells.y; ++j)
    for (int k = 0; k < numCells.z; ++k)
    {
        int idx = numCells.y * numCells.z * i + numCells.z * j + k;
        Int3 nodeIdx = *minCellIdx + Int3(i, j, k);
        b[idx].x = Bx(nodeIdx);
        b[idx].y = By(nodeIdx);
        b[idx].z = Bz(nodeIdx);
    }
}


void Grid::dumpE(FP3 * e, const Int3 * minCellIdx, const Int3 * maxCellIdx)
    const
{
    Int3 numCells = *maxCellIdx - *minCellIdx;
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < numCells.x; ++i)
    for (int j = 0; j < numCells.y; ++j)
    for (int k = 0; k < numCells.z; ++k)
    {
        int idx = numCells.y * numCells.z * i + numCells.z * j + k;
        Int3 nodeIdx = *minCellIdx + Int3(i, j, k);
        e[idx].x = Ex(nodeIdx);
        e[idx].y = Ey(nodeIdx);
        e[idx].z = Ez(nodeIdx);
    }
}

void Grid::dumpA(FP3 * a, const Int3 * minCellIdx,
        const Int3 * maxCellIdx) const
{
    Int3 numCells = *maxCellIdx - *minCellIdx;
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < numCells.x; ++i)
    for (int j = 0; j < numCells.y; ++j)
    for (int k = 0; k < numCells.z; ++k)
    {
        int idx = numCells.y * numCells.z * i + numCells.z * j + k;
        Int3 nodeIdx = *minCellIdx + Int3(i, j, k);
        a[idx].x = Ax(nodeIdx);
        a[idx].y = Ay(nodeIdx);
        a[idx].z = Az(nodeIdx);
    }
}

void Grid::dumpCurrents(FP3 * currents, const Int3 * minCellIdx,
    const Int3 * maxCellIdx) const
{
    Int3 numCells = *maxCellIdx - *minCellIdx;
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < numCells.x; ++i)
    for (int j = 0; j < numCells.y; ++j)
    for (int k = 0; k < numCells.z; ++k)
    {
        int idx = numCells.y * numCells.z * i + numCells.z * j + k;
        Int3 nodeIdx = *minCellIdx + Int3(i, j, k);
        currents[idx].x = Jx(nodeIdx);
        currents[idx].y = Jy(nodeIdx);
        currents[idx].z = Jz(nodeIdx);
        idx++;
    }
}

void Grid::loadE(const FP3 * e, const Int3 * minCellIdx, const Int3 * maxCellIdx)
{
    Int3 numCells = *maxCellIdx - *minCellIdx;
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < numCells.x; i++)
    for (int j = 0; j < numCells.y; j++)
    for (int k = 0; k < numCells.z; k++)
    {
        int idx = numCells.y * numCells.z * i + numCells.z * j + k;
        Int3 nodeIdx = *minCellIdx + Int3(i, j, k);
        Ex(nodeIdx) = e[idx].x;
        Ey(nodeIdx) = e[idx].y;
        Ez(nodeIdx) = e[idx].z;
    }
}

void Grid::loadB(const FP3 * b, const Int3 * minCellIdx, const Int3 * maxCellIdx)
{
    Int3 numCells = *maxCellIdx - *minCellIdx;
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < numCells.x; ++i)
    for (int j = 0; j < numCells.y; ++j)
    for (int k = 0; k < numCells.z; ++k)
    {
        int idx = numCells.y * numCells.z * i + numCells.z * j + k;
        Int3 nodeIdx = *minCellIdx + Int3(i, j, k);
        Bx(nodeIdx) = b[idx].x;
        By(nodeIdx) = b[idx].y;
        Bz(nodeIdx) = b[idx].z;
    }
}

void Grid::loadA(const FP3 * a, const Int3 * minCellIdx, const Int3 * maxCellIdx)
{
    Int3 numCells = *maxCellIdx - *minCellIdx;
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < numCells.x; i++)
    for (int j = 0; j < numCells.y; j++)
    for (int k = 0; k < numCells.z; k++)
    {
        int idx = numCells.y * numCells.z * i + numCells.z * j + k;
        Int3 nodeIdx = *minCellIdx + Int3(i, j, k);
        Ax(nodeIdx) = a[idx].x;
        Ay(nodeIdx) = a[idx].y;
        Az(nodeIdx) = a[idx].z;
    }
}

void Grid::loadCurrents(const FP3 * currents, const Int3 * minCellIdx, const Int3 * maxCellIdx)
{
    Int3 numCells = *maxCellIdx - *minCellIdx;
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < numCells.x; i++)
    for (int j = 0; j < numCells.y; j++)
    for (int k = 0; k < numCells.z; k++)
    {
        int idx = numCells.y * numCells.z * i + numCells.z * j + k;
        Int3 nodeIdx = *minCellIdx + Int3(i, j, k);
        Jx(nodeIdx) = currents[idx].x;
        Jy(nodeIdx) = currents[idx].y;
        Jz(nodeIdx) = currents[idx].z;
    }
}

} // namespace pica
