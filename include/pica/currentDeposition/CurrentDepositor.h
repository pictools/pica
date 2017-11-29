#ifndef PICA_CURRENT_DEPOSITOR_H
#define PICA_CURRENT_DEPOSITOR_H


#include "pica/grid/YeeGrid.h"


namespace pica {


template<class Grid>
class CurrentDepositorCIC {
public:
    CurrentDepositorCIC(Grid& grid);
    void deposit(const typename Grid::PositionType& position, const Vector3<typename Grid::ValueType>& current);
};



template<Dimension dimension, typename Real>
class CurrentDepositorCIC<YeeGrid<dimension, Real> > {
public:
    typedef YeeGrid<dimension, Real> GridType;
    typedef typename GridType::PositionType PositionType;
    typedef typename GridType::IndexType IndexType;
    typedef typename GridType::ValueType ValueType;
    typedef typename ScalarType<PositionType>::Type ScalarPositionType;

    CurrentDepositorCIC(GridType& grid) :
        grid(grid),
        inverseStep(inverse(grid.getStep())),
        normalizedOrigin(grid.getOrigin() / grid.getStep()),
        normalizedOriginStaggered(normalizedOrigin + ones<dimension, ScalarPositionType>() * static_cast<ScalarPositionType>(0.5)) {}

    void deposit(const PositionType& position, const Vector3<ValueType>& current)
    {
        IndexType indexCollocated, indexStaggered;
        PositionType coeffCollocated, coeffStaggered;
        getIndexCoeff(position, indexCollocated, coeffCollocated, indexStaggered, coeffStaggered);
        deposit(&GridType::jx, indexCollocated, coeffCollocated, current.x);
        deposit(&GridType::jy, indexCollocated, coeffCollocated, current.y);
        deposit(&GridType::jz, indexCollocated, coeffCollocated, current.z);
    }

private:

    typedef ValueType&(GridType::*FieldComponent1d)(int);
    void deposit(FieldComponent1d component, IndexType baseIndex, PositionType coeff, Real value) const
    {
        (grid.*component)(baseIndex.x) = (1.0 - coeff.x) * value;
        (grid.*component)(baseIndex.x + 1) = coeff.x * value;
    }

    typedef ValueType&(GridType::*FieldComponent2d)(int, int);
    void deposit(FieldComponent2d component, IndexType baseIndex, PositionType coeff, Real value) const
    {
        (grid.*component)(baseIndex.x, baseIndex.y) = (1.0 - coeff.x) * (1.0 - coeff.y) * value;
        (grid.*component)(baseIndex.x, baseIndex.y + 1) = (1.0 - coeff.x) * coeff.y * value;
        (grid.*component)(baseIndex.x + 1, baseIndex.y) = coeff.x * (1.0 - coeff.y) * value;
        (grid.*component)(baseIndex.x + 1, baseIndex.y + 1) = coeff.x * coeff.y * value;
    }

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

    GridType& grid;
    PositionType inverseStep;
    PositionType normalizedOrigin, normalizedOriginStaggered;

};


} // namespace pica


#endif
