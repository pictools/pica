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
    using typename Grid<dimension, Real>::ArrayType;

    YeeGrid(const PositionType& origin, const PositionType& step, const IndexType& size);

    ArrayType& ex();
    const ArrayType& ex() const;
    ArrayType& ey();
    const ArrayType& ey() const;
    ArrayType& ez();
    const ArrayType& ez() const;

    ArrayType& bx();
    const ArrayType& bx() const;
    ArrayType& by();
    const ArrayType& by() const;
    ArrayType& bz();
    const ArrayType& bz() const;

    ArrayType& jx();
    const ArrayType& jx() const;
    ArrayType& jy();
    const ArrayType& jy() const;
    ArrayType& jz();
    const ArrayType& jz() const;

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
    using typename Grid<One, Real>::ArrayType;

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

    ArrayType& ex() { return Grid<One, Real>::ex(); }
    const ArrayType& ex() const { return Grid<One, Real>::ex(); }
    ArrayType& ey() { return Grid<One, Real>::ey(); }
    const ArrayType& ey() const { return Grid<One, Real>::ey(); }
    ArrayType& ez() { return Grid<One, Real>::ez(); }
    const ArrayType& ez() const { return Grid<One, Real>::ez(); }

    ArrayType& bx() { return Grid<One, Real>::bx(); }
    const ArrayType& bx() const { return Grid<One, Real>::bx(); }
    ArrayType& by() { return Grid<One, Real>::by(); }
    const ArrayType& by() const { return Grid<One, Real>::by(); }
    ArrayType& bz() { return Grid<One, Real>::bz(); }
    const ArrayType& bz() const { return Grid<One, Real>::bz(); }

    ArrayType& jx() { return Grid<One, Real>::jx(); }
    const ArrayType& jx() const { return Grid<One, Real>::jx(); }
    ArrayType& jy() { return Grid<One, Real>::jy(); }
    const ArrayType& jy() const { return Grid<One, Real>::jy(); }
    ArrayType& jz() { return Grid<One, Real>::jz(); }
    const ArrayType& jz() const { return Grid<One, Real>::jz(); }

    // Get shifts of components from a cell origin
    PositionType getShiftEx() const { return static_cast<ValueType>(0); }
    PositionType getShiftEy() const { return static_cast<ValueType>(0.5 * this->step); }
    PositionType getShiftEz() const { return static_cast<ValueType>(0.5 * this->step); }
    PositionType getShiftBx() const { return static_cast<ValueType>(0.5 * this->step); }
    PositionType getShiftBy() const { return static_cast<ValueType>(0); }
    PositionType getShiftBz() const { return static_cast<ValueType>(0); }
    PositionType getShiftJx() const { return static_cast<ValueType>(0); }
    PositionType getShiftJy() const { return static_cast<ValueType>(0.5 * this->step); }
    PositionType getShiftJz() const { return static_cast<ValueType>(0.5 * this->step); }
};

// This specialization adds overloads for accessing by two ints in addition to standard Vector2<int>
template<typename Real>
class YeeGrid<Two, Real> : public Grid<Two, Real> {
public:

    using typename Grid<Two, Real>::IndexType;
    using typename Grid<Two, Real>::PositionType;
    using typename Grid<Two, Real>::ValueType;
    using typename Grid<Two, Real>::ArrayType;

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

    ArrayType& ex() { return Grid<Two, Real>::ex(); }
    const ArrayType& ex() const { return Grid<Two, Real>::ex(); }
    ArrayType& ey() { return Grid<Two, Real>::ey(); }
    const ArrayType& ey() const { return Grid<Two, Real>::ey(); }
    ArrayType& ez() { return Grid<Two, Real>::ez(); }
    const ArrayType& ez() const { return Grid<Two, Real>::ez(); }

    ArrayType& bx() { return Grid<Two, Real>::bx(); }
    const ArrayType& bx() const { return Grid<Two, Real>::bx(); }
    ArrayType& by() { return Grid<Two, Real>::by(); }
    const ArrayType& by() const { return Grid<Two, Real>::by(); }
    ArrayType& bz() { return Grid<Two, Real>::bz(); }
    const ArrayType& bz() const { return Grid<Two, Real>::bz(); }

    ArrayType& jx() { return Grid<Two, Real>::jx(); }
    const ArrayType& jx() const { return Grid<Two, Real>::jx(); }
    ArrayType& jy() { return Grid<Two, Real>::jy(); }
    const ArrayType& jy() const { return Grid<Two, Real>::jy(); }
    ArrayType& jz() { return Grid<Two, Real>::jz(); }
    const ArrayType& jz() const { return Grid<Two, Real>::jz(); }

    // Get shifts of components from a cell origin
    PositionType getShiftEx() const { return PositionType(0, 0.5 * this->step.y); }
    PositionType getShiftEy() const { return PositionType(0.5 * this->step.x, 0); }
    PositionType getShiftEz() const { return PositionType(0.5 * this->step.x, 0.5 * this->step.y); }
    PositionType getShiftBx() const { return PositionType(0.5 * this->step.x, 0); }
    PositionType getShiftBy() const { return PositionType(0, 0.5 * this->step.y); }
    PositionType getShiftBz() const { return PositionType(0, 0); }
    PositionType getShiftJx() const { return PositionType(0, 0.5 * this->step.y); }
    PositionType getShiftJy() const { return PositionType(0.5 * this->step.x, 0); }
    PositionType getShiftJz() const { return PositionType(0.5 * this->step.x, 0.5 * this->step.y); }
};

// This specialization adds overloads for accessing by two ints in addition to standard Vector3<int>
template<typename Real>
class YeeGrid<Three, Real> : public Grid<Three, Real> {
public:

    using typename Grid<Three, Real>::IndexType;
    using typename Grid<Three, Real>::PositionType;
    using typename Grid<Three, Real>::ValueType;
    using typename Grid<Three, Real>::ArrayType;

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

    ArrayType& ex() { return Grid<Three, Real>::ex(); }
    const ArrayType& ex() const { return Grid<Three, Real>::ex(); }
    ArrayType& ey() { return Grid<Three, Real>::ey(); }
    const ArrayType& ey() const { return Grid<Three, Real>::ey(); }
    ArrayType& ez() { return Grid<Three, Real>::ez(); }
    const ArrayType& ez() const { return Grid<Three, Real>::ez(); }

    ArrayType& bx() { return Grid<Three, Real>::bx(); }
    const ArrayType& bx() const { return Grid<Three, Real>::bx(); }
    ArrayType& by() { return Grid<Three, Real>::by(); }
    const ArrayType& by() const { return Grid<Three, Real>::by(); }
    ArrayType& bz() { return Grid<Three, Real>::bz(); }
    const ArrayType& bz() const { return Grid<Three, Real>::bz(); }

    ArrayType& jx() { return Grid<Three, Real>::jx(); }
    const ArrayType& jx() const { return Grid<Three, Real>::jx(); }
    ArrayType& jy() { return Grid<Three, Real>::jy(); }
    const ArrayType& jy() const { return Grid<Three, Real>::jy(); }
    ArrayType& jz() { return Grid<Three, Real>::jz(); }
    const ArrayType& jz() const { return Grid<Three, Real>::jz(); }

    // Get shifts of components from a cell origin
    PositionType getShiftEx() const { return PositionType(0, 0.5 * this->step.y, 0.5 * this->step.z); }
    PositionType getShiftEy() const { return PositionType(0.5 * this->step.x, 0, 0.5 * this->step.z); }
    PositionType getShiftEz() const { return PositionType(0.5 * this->step.x, 0.5 * this->step.y, 0); }
    PositionType getShiftBx() const { return PositionType(0.5 * this->step.x, 0, 0); }
    PositionType getShiftBy() const { return PositionType(0, 0.5 * this->step.y, 0); }
    PositionType getShiftBz() const { return PositionType(0, 0, 0.5 * this->step.z); }
    PositionType getShiftJx() const { return PositionType(0, 0.5 * this->step.y, 0.5 * this->step.z); }
    PositionType getShiftJy() const { return PositionType(0.5 * this->step.x, 0, 0.5 * this->step.z); }
    PositionType getShiftJz() const { return PositionType(0.5 * this->step.x, 0.5 * this->step.y, 0); }
};

} // namespace pica


#endif
