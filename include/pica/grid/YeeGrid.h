#ifndef PICA_YEE_GRID_H
#define PICA_YEE_GRID_H


#include "pica/grid/Grid.h"


namespace pica {


// In 3D index(i, j, k) refers to cell(i, j, k) being the cell between nodes(i, j, k)
//        and (i + 1, j + 1, k + 1).It corresponds to Ex(i, j + 1 / 2, k + 1 / 2),
//        Ey(i + 1 / 2, j, k + 1 / 2), Ez(i + 1 / 2, j + 1 / 2, k), same for the components of J,
//        Bx(i + 1 / 2, j, k), By(i, j + 1 / 2, k), Bz(i, j, k + 1 / 2).
template<Dimension dimension, typename Real = double>
class YeeGrid : public Grid<dimension, Real> {
public:

    using typename Grid<dimension, Real>::IndexType;
    using typename Grid<dimension, Real>::PositionType;

    YeeGrid(const PositionType& origin, const PositionType& step, const IndexType& size);

    PositionType getShiftEx() const;
    PositionType getShiftEy() const;
    PositionType getShiftEz() const;
    PositionType getShiftBx() const;
    PositionType getShiftBy() const;
    PositionType getShiftBz() const;
    PositionType getShiftJx() const;
    PositionType getShiftJy() const;
    PositionType getShiftJz() const;
};

// This specialization adds getShift-type functions
template<typename Real>
class YeeGrid<One, Real> : public Grid<One, Real> {
public:

    using typename Grid<One, Real>::IndexType;
    using typename Grid<One, Real>::PositionType;
    using typename Grid<One, Real>::ValueType;

    YeeGrid(const PositionType& origin, const PositionType& step, const IndexType& size) :
        Grid<One, Real>(origin, step, size) {}

    ValueType& ex(int i) { return Grid<One, Real>::ex(IndexType(i)); }
    ValueType ex(int i) const { return Grid<One, Real>::ex(IndexType(i)); }
    ValueType& ey(int i) { return Grid<One, Real>::ey(IndexType(i)); }
    ValueType ey(int i) const { return Grid<One, Real>::ey(IndexType(i)); }
    ValueType& ez(int i) { return Grid<One, Real>::ez(IndexType(i)); }
    ValueType ez(int i) const { return Grid<One, Real>::ez(IndexType(i)); }

    ValueType& bx(int i) { return Grid<One, Real>::bx(IndexType(i)); }
    ValueType bx(int i) const { return Grid<One, Real>::bx(IndexType(i)); }
    ValueType& by(int i) { return Grid<One, Real>::by(IndexType(i)); }
    ValueType by(int i) const { return Grid<One, Real>::by(IndexType(i)); }
    ValueType& bz(int i) { return Grid<One, Real>::bz(IndexType(i)); }
    ValueType bz(int i) const { return Grid<One, Real>::bz(IndexType(i)); }

    ValueType& jx(int i) { return Grid<One, Real>::jx(IndexType(i)); }
    ValueType jx(int i) const { return Grid<One, Real>::jx(IndexType(i)); }
    ValueType& jy(int i) { return Grid<One, Real>::jy(IndexType(i)); }
    ValueType jy(int i) const { return Grid<One, Real>::jy(IndexType(i)); }
    ValueType& jz(int i) { return Grid<One, Real>::jz(IndexType(i)); }
    ValueType jz(int i) const { return Grid<One, Real>::jz(IndexType(i)); }

    // Get shifts of components from a cell origin
    PositionType getShiftEx() const { return static_cast<ValueType>(0); }
    PositionType getShiftEy() const { return static_cast<ValueType>(0.5 * step); }
    PositionType getShiftEz() const { return static_cast<ValueType>(0.5 * step); }
    PositionType getShiftBx() const { return static_cast<ValueType>(0.5 * step); }
    PositionType getShiftBy() const { return static_cast<ValueType>(0); }
    PositionType getShiftBz() const { return static_cast<ValueType>(0); }
    PositionType getShiftJx() const { return static_cast<ValueType>(0); }
    PositionType getShiftJy() const { return static_cast<ValueType>(0.5 * step); }
    PositionType getShiftJz() const { return static_cast<ValueType>(0.5 * step); }
};

// This specialization adds overloads for accessing by two ints in addition to standard Vector2<int>
template<typename Real>
class YeeGrid<Two, Real> : public Grid<Two, Real> {
public:

    using typename Grid<Two, Real>::IndexType;
    using typename Grid<Two, Real>::PositionType;
    using typename Grid<Two, Real>::ValueType;

    YeeGrid(const PositionType& origin, const PositionType& step, const IndexType& size) :
        Grid<Two, Real>(origin, step, size) {}

    ValueType& ex(int i, int j) { return Grid<Two, Real>::ex(IndexType(i, j)); }
    ValueType ex(int i, int j) const { return Grid<Two, Real>::ex(IndexType(i, j)); }
    ValueType& ey(int i, int j) { return Grid<Two, Real>::ey(IndexType(i, j)); }
    ValueType ey(int i, int j) const { return Grid<Two, Real>::ey(IndexType(i, j)); }
    ValueType& ez(int i, int j) { return Grid<Two, Real>::ez(IndexType(i, j)); }
    ValueType ez(int i, int j) const { return Grid<Two, Real>::ez(IndexType(i, j)); }

    ValueType& bx(int i, int j) { return Grid<Two, Real>::bx(IndexType(i, j)); }
    ValueType bx(int i, int j) const { return Grid<Two, Real>::bx(IndexType(i, j)); }
    ValueType& by(int i, int j) { return Grid<Two, Real>::by(IndexType(i, j)); }
    ValueType by(int i, int j) const { return Grid<Two, Real>::by(IndexType(i, j)); }
    ValueType& bz(int i, int j) { return Grid<Two, Real>::bz(IndexType(i, j)); }
    ValueType bz(int i, int j) const { return Grid<Two, Real>::bz(IndexType(i, j)); }

    ValueType& jx(int i, int j) { return Grid<Two, Real>::jx(IndexType(i, j)); }
    ValueType jx(int i, int j) const { return Grid<Two, Real>::jx(IndexType(i, j)); }
    ValueType& jy(int i, int j) { return Grid<Two, Real>::jy(IndexType(i, j)); }
    ValueType jy(int i, int j) const { return Grid<Two, Real>::jy(IndexType(i, j)); }
    ValueType& jz(int i, int j) { return Grid<Two, Real>::jz(IndexType(i, j)); }
    ValueType jz(int i, int j) const { return Grid<Two, Real>::jz(IndexType(i, j)); }

    // Get shifts of components from a cell origin
    PositionType getShiftEx() const { return PositionType(0, 0.5 * step.y); }
    PositionType getShiftEy() const { return PositionType(0.5 * step.x, 0); }
    PositionType getShiftEz() const { return PositionType(0.5 * step.x, 0.5 * step.y); }
    PositionType getShiftBx() const { return PositionType(0.5 * step.x, 0); }
    PositionType getShiftBy() const { return PositionType(0, 0.5 * step.y); }
    PositionType getShiftBz() const { return PositionType(0, 0); }
    PositionType getShiftJx() const { return PositionType(0, 0.5 * step.y); }
    PositionType getShiftJy() const { return PositionType(0.5 * step.x, 0); }
    PositionType getShiftJz() const { return PositionType(0.5 * step.x, 0.5 * step.y); }
};

// This specialization adds overloads for accessing by two ints in addition to standard Vector3<int>
template<typename Real>
class YeeGrid<Three, Real> : public Grid<Three, Real> {
public:

    using typename Grid<Three, Real>::IndexType;
    using typename Grid<Three, Real>::PositionType;
    using typename Grid<Three, Real>::ValueType;

    YeeGrid(const PositionType& origin, const PositionType& step, const IndexType& size) :
        Grid<Three, Real>(origin, step, size) {}

    ValueType& ex(int i, int j, int k) { return Grid<Three, Real>::ex(IndexType(i, j, k)); }
    ValueType ex(int i, int j, int k) const { return Grid<Three, Real>::ex(IndexType(i, j, k)); }
    ValueType& ey(int i, int j, int k) { return Grid<Three, Real>::ey(IndexType(i, j, k)); }
    ValueType ey(int i, int j, int k) const { return Grid<Three, Real>::ey(IndexType(i, j, k)); }
    ValueType& ez(int i, int j, int k) { return Grid<Three, Real>::ez(IndexType(i, j, k)); }
    ValueType ez(int i, int j, int k) const { return Grid<Three, Real>::ez(IndexType(i, j, k)); }

    ValueType& bx(int i, int j, int k) { return Grid<Three, Real>::bx(IndexType(i, j, k)); }
    ValueType bx(int i, int j, int k) const { return Grid<Three, Real>::bx(IndexType(i, j, k)); }
    ValueType& by(int i, int j, int k) { return Grid<Three, Real>::by(IndexType(i, j, k)); }
    ValueType by(int i, int j, int k) const { return Grid<Three, Real>::by(IndexType(i, j, k)); }
    ValueType& bz(int i, int j, int k) { return Grid<Three, Real>::bz(IndexType(i, j, k)); }
    ValueType bz(int i, int j, int k) const { return Grid<Three, Real>::bz(IndexType(i, j, k)); }

    ValueType& jx(int i, int j, int k) { return Grid<Three, Real>::jx(IndexType(i, j, k)); }
    ValueType jx(int i, int j, int k) const { return Grid<Three, Real>::jx(IndexType(i, j, k)); }
    ValueType& jy(int i, int j, int k) { return Grid<Three, Real>::jy(IndexType(i, j, k)); }
    ValueType jy(int i, int j, int k) const { return Grid<Three, Real>::jy(IndexType(i, j, k)); }
    ValueType& jz(int i, int j, int k) { return Grid<Three, Real>::jz(IndexType(i, j, k)); }
    ValueType jz(int i, int j, int k) const { return Grid<Three, Real>::jz(IndexType(i, j, k)); }

    // Get shifts of components from a cell origin
    PositionType getShiftEx() const { return PositionType(0, 0.5 * step.y, 0.5 * step.z); }
    PositionType getShiftEy() const { return PositionType(0.5 * step.x, 0, 0.5 * step.z); }
    PositionType getShiftEz() const { return PositionType(0.5 * step.x, 0.5 * step.y, 0); }
    PositionType getShiftBx() const { return PositionType(0.5 * step.x, 0, 0); }
    PositionType getShiftBy() const { return PositionType(0, 0.5 * step.y, 0); }
    PositionType getShiftBz() const { return PositionType(0, 0, 0.5 * step.z); }
    PositionType getShiftJx() const { return PositionType(0, 0.5 * step.y, 0.5 * step.z); }
    PositionType getShiftJy() const { return PositionType(0.5 * step.x, 0, 0.5 * step.z); }
    PositionType getShiftJz() const { return PositionType(0.5 * step.x, 0.5 * step.y, 0); }
};

} // namespace pica


#endif
