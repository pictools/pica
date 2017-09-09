#ifndef PICA_GRID_H
#define PICA_GRID_H


#include "pica/grid/ScalarField.h"
#include "pica/math/Vectors.h"
#include "pica/utility/Array.h"


namespace pica {


template<Dimension dimension, typename Real = double>
class Grid_ {
public:
    typedef Real ValueType;
    typedef typename VectorTypeHelper<dimension, Real>::Type PositionType;
    typedef typename VectorTypeHelper<dimension, int>::Type IndexType;

    Grid_(const PositionType& origin, const PositionType& step, const IndexType& size) :
        exData(size),
        eyData(size),
        ezData(size),
        bxData(size),
        byData(size),
        bzData(size),
        jxData(size),
        jyData(size),
        jzData(size),
        origin(origin),
        step(step) {}

    IndexType getSize() const { return exData.getSize(); }
    PositionType getOrigin() const { return origin; }
    PositionType getStep() const { return step; }
    IndexType getCellIndex(const PositionType& position) { return truncate((position - origin) / step); }

    ValueType& ex(const IndexType& index) { return exData(index); }
    ValueType ex(const IndexType& index) const { return exData(index); }
    ValueType& ey(const IndexType& index) { return eyData(index); }
    ValueType ey(const IndexType& index) const { return eyData(index); }
    ValueType& ez(const IndexType& index) { return ezData(index); }
    ValueType ez(const IndexType& index) const { return ezData(index); }

    ValueType& bx(const IndexType& index) { return bxData(index); }
    ValueType bx(const IndexType& index) const { return bxData(index); }
    ValueType& by(const IndexType& index) { return byData(index); }
    ValueType by(const IndexType& index) const { return byData(index); }
    ValueType& bz(const IndexType& index) { return bzData(index); }
    ValueType bz(const IndexType& index) const { return bzData(index); }

    ValueType& jx(const IndexType& index) { return jxData(index); }
    ValueType jx(const IndexType& index) const { return jxData(index); }
    ValueType& jy(const IndexType& index) { return jyData(index); }
    ValueType jy(const IndexType& index) const { return jyData(index); }
    ValueType& jz(const IndexType& index) { return jzData(index); }
    ValueType jz(const IndexType& index) const { return jzData(index); }

private:
    typedef typename ArrayTypeHelper<dimension, Real>::Type ArrayType;
    ArrayType exData, eyData, ezData, bxData, byData, bzData, jxData, jyData, jzData;
    PositionType origin, step;
};



/* Yee grid for FDTD method. Implementation is based on Computational
Electrodynamics: The Finite-Difference Time-Domain Method by Allen Taflove.
The class provides index-based access and trilinear interpolation from 8
nearest grid values.
Components of the electric and magnetic fields correspond to different points.
Index (i, j, k) refers to cell (i, j, k) being the cell between nodes (i, j, k)
and (i + 1, j + 1, k + 1). It corresponds to Ex(i, j+1/2, k+1/2),
Ey(i+1/2, j, k+1/2), Ez(i+1/2, j+1/2, k), same for the components of J,
Bx(i+1/2, j, k), By(i, j+1/2, k), Bz(i, j, k+1/2).
There is also EB_timeShift = dt/2 time difference between E(J) and B as needed
for FDTD.
For field exchanges and interpolation we need to store some grid values
that correspond to points outside of the main grid area, size of internal and
external area are given as constant members. */
class Grid
{

public:

    Grid(const Int3 & _numInternalCells, FP _dt,
        const FP3 & minCoords, const FP3 & _steps,
        const Int3 & globalGridDims, const Int3 & domainOffset);

    const FP3 BxPosition(int x, int y, int z) const
    {
        return baseCoords(x, y, z) + shiftBx;
    }

    const FP3 ByPosition(int x, int y, int z) const
    {
        return baseCoords(x, y, z) + shiftBy;
    }

    const FP3 BzPosition(int x, int y, int z) const
    {
        return baseCoords(x, y, z) + shiftBz;
    }

    const FP3 ExPosition(int x, int y, int z) const
    {
        return baseCoords(x, y, z) + shiftEJx;
    }

    const FP3 EyPosition(int x, int y, int z) const
    {
        return baseCoords(x, y, z) + shiftEJy;
    }

    const FP3 EzPosition(int x, int y, int z) const
    {
        return baseCoords(x, y, z) + shiftEJz;
    }

    const FP3 JxPosition(int x, int y, int z) const
    {
        return baseCoords(x, y, z) + shiftEJx;
    }

    const FP3 JyPosition(int x, int y, int z) const
    {
        return baseCoords(x, y, z) + shiftEJy;
    }

    const FP3 JzPosition(int x, int y, int z) const
    {
        return baseCoords(x, y, z) + shiftEJz;
    }

    const FP3 RhoPosition(int x, int y, int z) const
    {
        return baseCoords(x, y, z) + shiftRho;
    }

    void getField(FP x, FP y, FP z, FP3 & e, FP3 & b) const
    {
        getFields(x, y, z, e, b);
    }
    void getFields(FP x, FP y, FP z, FP3 & e, FP3 & b) const;
    void getJ(FP x, FP y, FP z, FP3 & j) const;
    void getA(FP x, FP y, FP z, FP3 & a) const;

    void getFieldsCIC(FP x, FP y, FP z, FP3 & e, FP3 & b) const;
    void getFieldsTSC(FP x, FP y, FP z, FP3 & e, FP3 & b) const;
    void getFieldsSecondOrder(FP x, FP y, FP z, FP3 & e, FP3 & b) const;
    void getFieldsFourthOrder(FP x, FP y, FP z, FP3 & e, FP3 & b) const;
    void getFieldsPCS(FP x, FP y, FP z, FP3 & e, FP3 & b) const;

    // Single-component interpolation for graphics
    FP getExCIC(const FP3& coords) const;
    FP getEyCIC(const FP3& coords) const;
    FP getEzCIC(const FP3& coords) const;
    FP getBxCIC(const FP3& coords) const;
    FP getByCIC(const FP3& coords) const;
    FP getBzCIC(const FP3& coords) const;
    FP getJxCIC(const FP3& coords) const;
    FP getJyCIC(const FP3& coords) const;
    FP getJzCIC(const FP3& coords) const;
    FP getAxCIC(const FP3& coords) const;
    FP getAyCIC(const FP3& coords) const;
    FP getAzCIC(const FP3& coords) const;
    FP getExSecondOrder(const FP3& coords) const;
    FP getEySecondOrder(const FP3& coords) const;
    FP getEzSecondOrder(const FP3& coords) const;
    FP getBxSecondOrder(const FP3& coords) const;
    FP getBySecondOrder(const FP3& coords) const;
    FP getBzSecondOrder(const FP3& coords) const;
    FP getJxSecondOrder(const FP3& coords) const;
    FP getJySecondOrder(const FP3& coords) const;
    FP getJzSecondOrder(const FP3& coords) const;
    FP getAxSecondOrder(const FP3& coords) const;
    FP getAySecondOrder(const FP3& coords) const;
    FP getAzSecondOrder(const FP3& coords) const;

    void depositJ(const FP3 & value, const FP3 & coords);
    void depositJCIC(const FP3 & value, const FP3 & coords);
    void depositJTSC(const FP3 & value, const FP3 & coords);
    void depositJPCS(const FP3 & value, const FP3 & coords);

    void addCurrents(const FP3 * currents, const Int3 * minCellIdx, const Int3 * maxCellIdx);
    void dumpE(FP3 * e, const Int3 * minCellIdx, const Int3 * maxCellIdx) const;
    void dumpB(FP3 * b, const Int3 * minCellIdx, const Int3 * maxCellIdx) const;
    void dumpA(FP3 * a, const Int3 * minCellIdx,  const Int3 * maxCellIdx) const;
    void dumpCurrents(FP3 * currents, const Int3 * minCellIdx, const Int3 * maxCellIdx) const;

    virtual void loadE(const FP3 * e, const Int3 * minCellIdx, const Int3 * maxCellIdx);
    virtual void loadB(const FP3 * b, const Int3 * minCellIdx, const Int3 * maxCellIdx);
    virtual void loadA(const FP3 * a, const Int3 * minCellIdx, const Int3 * maxCellIdx);
    virtual void loadCurrents(const FP3 * currents, const Int3 * minCellIdx, const Int3 * maxCellIdx);

    /* Make all current density values zero. */
    void zeroizeJ();

    const Int3 getNumExternalLeftCells() const
    {
        Int3 result(2, 2, 2);
        for (int d = 0; d < 3; d++)
            if (globalGridDims[d] == 1)
                result[d] = 0;
        return result;
    }

    const Int3 getNumExternalRightCells() const
    {
        return getNumExternalLeftCells();
    }

    enum InterpolationType { Interpolation_CIC, Interpolation_TSC,
        Interpolation_SecondOrder, Interpolation_FourthOrder, Interpolation_PCS };
    void setInterpolationType(InterpolationType type);
    InterpolationType getInterpolationType() const;

    enum DepositionType { Deposition_CIC, Deposition_VillasenorBuneman,
        Deposition_ZigzagFirstOrder, Deposition_TSC,
        Deposition_ZigzagSecondOrder, Deposition_Esirkepov, Deposition_PCS };
    void setDepositionType(DepositionType type);
    DepositionType getDepositionType() const;

    const Int3 globalGridDims, domainOffset; // important to initialize it first
    const FP3 steps;
    const FP dt;
    const Int3 numInternalCells;
    const Int3 numCells;
    const FP3 origin;

    // Time diffence between b and e
    const FP EB_timeShift;

    const FP3 shiftRho;

    ScalarField Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz;
    ScalarField Ax, Ay, Az; // vector potential

private:

    // 3d shifts of the field in the cell
    const FP3 shiftEJx, shiftEJy, shiftEJz,
        shiftBx, shiftBy, shiftBz;

    /* Get grid index and normalized internal coords in [0, 0, 0]..(1, 1, 1) for
       given physical coords and shift. */
    void getGridCoords(const FP3 & coords, const FP3 & shift, Int3 & idx,
        FP3 & internalCoords) const
    {
        idx.x = (int)((coords.x - origin.x - shift.x) / steps.x);
        idx.y = (int)((coords.y - origin.y - shift.y) / steps.y);
        idx.z = (int)((coords.z - origin.z - shift.z) / steps.z);
        internalCoords = (coords - baseCoords(idx.x, idx.y, idx.z) - shift) / steps;
    }

    void getClosestGridCoords(const FP3 & coords, const FP3 & shift, Int3 & idx,
        FP3 & internalCoords) const
    {
        idx.x = (int)((coords.x - origin.x - shift.x) / steps.x + 0.5);
        idx.y = (int)((coords.y - origin.y - shift.y) / steps.y + 0.5);
        idx.z = (int)((coords.z - origin.z - shift.z) / steps.z + 0.5);
        internalCoords = (coords - baseCoords(idx.x, idx.y, idx.z) - shift) / steps;
    }

    /* Get base coords of element (i, j, k) so that its real coords are
       base coords + corresponding shift. */
    const FP3 baseCoords(int i, int j, int k) const
    {
        return origin + FP3(i, j, k) * steps;
    }

    InterpolationType interpolationType;
    void (Grid::*interpolationFunction)(FP, FP, FP, FP3&, FP3&) const;
    DepositionType depositionType;
    void (Grid::*depositionFunction)(const FP3&, const FP3&);
};

} // namespace pica


#endif
