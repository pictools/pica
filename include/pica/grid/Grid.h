#ifndef PICA_GRID_H
#define PICA_GRID_H


#include "pica/math/Vectors.h"
#include "pica/utility/Array.h"


namespace pica {


template<Dimension dimension, typename Real = double>
class Grid {
public:
    typedef Real ValueType;
    typedef typename VectorTypeHelper<dimension, Real>::Type PositionType;
    typedef typename VectorTypeHelper<dimension, int>::Type IndexType;
    typedef typename ArrayTypeHelper<dimension, Real>::Type ArrayType;

    Grid(const PositionType& origin, const PositionType& step, const IndexType& size) :
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

    ArrayType& ex() { return exData; }
    const ArrayType& ex() const { return exData; }
    ArrayType& ey() { return eyData; }
    const ArrayType& ey() const { return eyData; }
    ArrayType& ez() { return ezData; }
    const ArrayType& ez() const { return ezData; }

    ArrayType& bx() { return bxData; }
    const ArrayType& bx() const { return bxData; }
    ArrayType& by() { return byData; }
    const ArrayType& by() const { return byData; }
    ArrayType& bz() { return bzData; }
    const ArrayType& bz() const { return bzData; }

    ArrayType& jx() { return jxData; }
    const ArrayType& jx() const { return jxData; }
    ArrayType& jy() { return jyData; }
    const ArrayType& jy() const { return jyData; }
    ArrayType& jz() { return jzData; }
    const ArrayType& jz() const { return jzData; }

private:
    ArrayType exData, eyData, ezData, bxData, byData, bzData, jxData, jyData, jzData;
    PositionType origin, step;
};


} // namespace pica


#endif
