#ifndef PICA_CURRENT_DEPOSITOR_YEE_GRID_H
#define PICA_CURRENT_DEPOSITOR_YEE_GRID_H


#include "pica/currentDeposition/CurrentDepositorYeeGrid.h"
#include "pica/grid/YeeGrid.h"


namespace pica {


namespace internal {

template<Dimension dimension, typename Real>
class CurrentDepositorYeeGridCICBase {
public:
    typedef YeeGrid<dimension, Real> GridType;
    typedef typename GridType::PositionType PositionType;
    typedef typename GridType::IndexType IndexType;
    typedef typename GridType::ValueType ValueType;
    typedef typename ScalarType<PositionType>::Type ScalarPositionType;

    CurrentDepositorYeeGridCICBase(GridType& grid) :
        grid(grid),
        inverseStep(inverse(grid.getStep())),
        normalizedOrigin(grid.getOrigin() / grid.getStep()),
        normalizedOriginStaggered(normalizedOrigin + ones<dimension, ScalarPositionType>() * static_cast<ScalarPositionType>(0.5)) {}

    void getIndexCoeff(const PositionType& position, IndexType& indexCollocated, PositionType& coeffCollocated,
        IndexType& indexStaggered, PositionType& coeffStaggered) const
    {
        PositionType normalizedPosition = position * inverseStep;
        PositionType internalPosition = normalizedPosition - normalizedOrigin;
        indexCollocated = truncate(internalPosition);
        coeffCollocated = internalPosition - PositionType(indexCollocated);
        internalPosition = normalizedPosition - normalizedOriginStaggered;
        indexStaggered = truncate(internalPosition);
        coeffStaggered = internalPosition - PositionType(indexStaggered);
    }

protected:
    GridType& grid;
    PositionType inverseStep;
    PositionType normalizedOrigin, normalizedOriginStaggered;

};    

} // namespace pica::internal


template<typename Real>
class CurrentDepositorCIC<YeeGrid<One, Real> > : public internal::CurrentDepositorYeeGridCICBase<One, Real> {
public:

    CurrentDepositorCIC(GridType& grid) : internal::CurrentDepositorYeeGridCICBase<One, Real>(grid)
    {}
 
    void deposit(const PositionType& position, const Vector3<ValueType>& current)
    {
        IndexType indexCollocated, indexStaggered;
        PositionType coeffCollocated, coeffStaggered;
        getIndexCoeff(position, indexCollocated, coeffCollocated, indexStaggered, coeffStaggered);
        deposit(&GridType::jx, IndexType(indexCollocated.x), PositionType(coeffCollocated.x), current.x);
        deposit(&GridType::jy, IndexType(indexStaggered.x), PositionType(coeffStaggered.x), current.y);
        deposit(&GridType::jz, IndexType(indexStaggered.x), PositionType(coeffStaggered.x), current.z);
    }

private:

    typedef ValueType&(GridType::*FieldComponent1d)(int);
    void deposit(FieldComponent1d component, IndexType baseIndex, PositionType coeff, Real value) const
    {
        (grid.*component)(baseIndex.x) = (1.0 - coeff.x) * value;
        (grid.*component)(baseIndex.x + 1) = coeff.x * value;
    }

};


template<typename Real>
class CurrentDepositorCIC<YeeGrid<Two, Real> > : public internal::CurrentDepositorYeeGridCICBase<Two, Real> {
public:

    CurrentDepositorCIC(GridType& grid) : internal::CurrentDepositorYeeGridCICBase<Two, Real>(grid)
    {}

    void deposit(const PositionType& position, const Vector3<ValueType>& current)
    {
        IndexType indexCollocated, indexStaggered;
        PositionType coeffCollocated, coeffStaggered;
        getIndexCoeff(position, indexCollocated, coeffCollocated, indexStaggered, coeffStaggered);
        deposit(&GridType::jx, IndexType(indexCollocated.x, indexStaggered.y),
            PositionType(coeffCollocated.x, coeffStaggered.y), current.x);
        deposit(&GridType::jy, IndexType(indexStaggered.x, indexCollocated.y),
            PositionType(coeffStaggered.x, coeffCollocated.y), current.y);
        deposit(&GridType::jz, IndexType(indexStaggered.x, indexStaggered.y),
            PositionType(coeffStaggered.x, coeffStaggered.y), current.z);
    }

private:

    typedef ValueType&(GridType::*FieldComponent2d)(int, int);
    void deposit(FieldComponent2d component, IndexType baseIndex, PositionType coeff, Real value) const
    {
        (grid.*component)(baseIndex.x, baseIndex.y) = (1.0 - coeff.x) * (1.0 - coeff.y) * value;
        (grid.*component)(baseIndex.x, baseIndex.y + 1) = (1.0 - coeff.x) * coeff.y * value;
        (grid.*component)(baseIndex.x + 1, baseIndex.y) = coeff.x * (1.0 - coeff.y) * value;
        (grid.*component)(baseIndex.x + 1, baseIndex.y + 1) = coeff.x * coeff.y * value;
    }

};


template<typename Real>
class CurrentDepositorCIC<YeeGrid<Three, Real> > : public internal::CurrentDepositorYeeGridCICBase<Three, Real> {
public:

    CurrentDepositorCIC(GridType& grid) : internal::CurrentDepositorYeeGridCICBase<Three, Real>(grid)
    {}

    void deposit(const PositionType& position, const Vector3<ValueType>& current)
    {
        IndexType indexCollocated, indexStaggered;
        PositionType coeffCollocated, coeffStaggered;
        getIndexCoeff(position, indexCollocated, coeffCollocated, indexStaggered, coeffStaggered);
        deposit(&GridType::jx, IndexType(indexCollocated.x, indexStaggered.y, indexStaggered.z),
            PositionType(coeffCollocated.x, coeffStaggered.y, coeffStaggered.z), current.x);
        deposit(&GridType::jy, IndexType(indexStaggered.x, indexCollocated.y, indexStaggered.z),
            PositionType(coeffStaggered.x, coeffCollocated.y, coeffStaggered.z), current.y);
        deposit(&GridType::jz, IndexType(indexStaggered.x, indexStaggered.y, indexCollocated.z),
            PositionType(coeffStaggered.x, coeffStaggered.y, coeffCollocated.z), current.z);
    }

private:

    typedef ValueType&(GridType::*FieldComponent3d)(int, int, int);
    void deposit(FieldComponent3d component, IndexType baseIndex, PositionType coeff, Real value) const
    {
        (grid.*component)(baseIndex.x, baseIndex.y, baseIndex.z) = (1.0 - coeff.x) * (1.0 - coeff.y) * (1.0 - coeff.z) * value;
        (grid.*component)(baseIndex.x, baseIndex.y, baseIndex.z + 1) = (1.0 - coeff.x) * (1.0 - coeff.y) * coeff.z * value;
        (grid.*component)(baseIndex.x, baseIndex.y + 1, baseIndex.z) = (1.0 - coeff.x) * coeff.y * (1.0 - coeff.z) * value;
        (grid.*component)(baseIndex.x, baseIndex.y + 1, baseIndex.z + 1) = (1.0 - coeff.x) * coeff.y * coeff.z * value;
        (grid.*component)(baseIndex.x + 1, baseIndex.y, baseIndex.z) = (1.0 - coeff.x) * (1.0 - coeff.y) * (1.0 - coeff.z) * value;
        (grid.*component)(baseIndex.x + 1, baseIndex.y, baseIndex.z + 1) = coeff.x * (1.0 - coeff.y) * coeff.z * value;
        (grid.*component)(baseIndex.x + 1, baseIndex.y + 1, baseIndex.z) = coeff.x * coeff.y * (1.0 - coeff.z) * value;
        (grid.*component)(baseIndex.x + 1, baseIndex.y + 1, baseIndex.z + 1) = coeff.x * coeff.y * coeff.z * value;;
    }

};


} // namespace pica


#endif
