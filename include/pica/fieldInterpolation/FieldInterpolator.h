#ifndef PICA_FIELD_INTERPOLATOR_H
#define PICA_FIELD_INTERPOLATOR_H


#include "pica/fieldInterpolation/Formfactor.h"
#include "pica/grid/YeeGrid.h"


namespace pica {


template<class Grid>
class FieldInterpolatorCIC {
public:
    FieldInterpolatorCIC(const Grid& grid);
    void get(const typename Grid::PositionType& position, Vector3<typename Grid::ValueType>& e, Vector3<typename Grid::ValueType>& b);
};



template<Dimension dimension, typename Real>
class FieldInterpolatorCIC<YeeGrid<dimension, Real> > {
public:
    typedef YeeGrid<dimension, Real> GridType;
    typedef typename GridType::PositionType PositionType;
    typedef typename GridType::IndexType IndexType;
    typedef typename GridType::ValueType ValueType;
    typedef typename ScalarType<PositionType>::Type ScalarPositionType;

    FieldInterpolatorCIC(const GridType& grid) :
        grid(grid),
        inverseStep(inverse(grid.getStep())),
        normalizedOrigin(grid.getOrigin() / grid.getStep()),
        normalizedOriginStaggered(normalizedOrigin + ones<dimension, ScalarPositionType>() * static_cast<ScalarPositionType>(0.5)) {}

    void get(const PositionType position, Vector3<ValueType>& e, Vector3<ValueType>& b)
    {
        IndexType indexCollocated, indexStaggered;
        PositionType coeffCollocated, coeffStaggered;
        getIndexCoeff(position, indexCollocated, coeffCollocated, indexStaggered, coeffStaggered);
        e.x = interpolate(&GridType::ex, indexCollocated, coeffCollocated);
        e.y = interpolate(&GridType::ey, indexStaggered, coeffStaggered);
        e.z = interpolate(&GridType::ez, indexStaggered, coeffStaggered);
        b.x = interpolate(&GridType::bx, indexStaggered, coeffStaggered);
        b.y = interpolate(&GridType::by, indexCollocated, coeffStaggered);
        b.z = interpolate(&GridType::bz, indexCollocated, coeffCollocated);
    }

private:

    typedef ValueType(GridType::*FieldComponent1d)(int) const;
    ValueType interpolate(FieldComponent1d component, IndexType baseIndex, PositionType coeff) const
    {
        return (grid.*component)(baseIndex.x) * (1.0 - coeff.x) +
            (grid.*component)(baseIndex.x + 1) * coeff.x;
    }

    typedef ValueType(GridType::*FieldComponent2d)(int, int) const;
    ValueType interpolate(FieldComponent2d component, const IndexType& baseIndex, const PositionType& coeff) const
    {
        return
            (grid.*component)(baseIndex.x, baseIndex.y) * (1.0 - coeff.x) * (1.0 - coeff.y) +
            (grid.*component)(baseIndex.x, baseIndex.y + 1) * (1.0 - coeff.x) * coeff.y +
            (grid.*component)(baseIndex.x + 1, baseIndex.y) * coeff.x * (1.0 - coeff.y) +
            (grid.*component)(baseIndex.x + 1, baseIndex.y + 1) * coeff.x * coeff.y;
    }

    typedef ValueType(GridType::*FieldComponent3d)(int, int, int) const;
    ValueType interpolate(FieldComponent3d component, const IndexType& baseIndex, const PositionType& coeff) const
    {
        return
            (grid.*component)(baseIndex.x, baseIndex.y, baseIndex.z) * (1.0 - coeff.x) * (1.0 - coeff.y) * (1.0 - coeff.z) +
            (grid.*component)(baseIndex.x, baseIndex.y, baseIndex.z + 1) * (1.0 - coeff.x) * (1.0 - coeff.y) * coeff.z +
            (grid.*component)(baseIndex.x, baseIndex.y + 1, baseIndex.z) * (1.0 - coeff.x) * coeff.y * (1.0 - coeff.z) +
            (grid.*component)(baseIndex.x, baseIndex.y + 1, baseIndex.z + 1) * (1.0 - coeff.x) * coeff.y * coeff.z +
            (grid.*component)(baseIndex.x + 1, baseIndex.y, baseIndex.z) * (1.0 - coeff.x) * (1.0 - coeff.y) * (1.0 - coeff.z) +
            (grid.*component)(baseIndex.x + 1, baseIndex.y, baseIndex.z + 1) * coeff.x * (1.0 - coeff.y) * coeff.z +
            (grid.*component)(baseIndex.x + 1, baseIndex.y + 1, baseIndex.z) * coeff.x * coeff.y * (1.0 - coeff.z) +
            (grid.*component)(baseIndex.x + 1, baseIndex.y + 1, baseIndex.z + 1) * coeff.x * coeff.y * coeff.z;
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

    const GridType& grid;
    PositionType inverseStep;
    PositionType normalizedOrigin, normalizedOriginStaggered;

};


} // namespace pica


#endif
