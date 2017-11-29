#ifndef PICA_FIELD_INTERPOLATOR_H
#define PICA_FIELD_INTERPOLATOR_H


#include "pica/grid/YeeGrid.h"


namespace pica {


template<class Grid>
class FieldInterpolatorCIC {
public:
    FieldInterpolatorCIC(const Grid& grid);
    void get(const typename Grid::PositionType& position, Vector3<typename Grid::ValueType>& e, Vector3<typename Grid::ValueType>& b);
};


namespace internal {

template<Dimension dimension, typename Real>
class FieldInterpolatorCICBase {
public:
    typedef YeeGrid<dimension, Real> GridType;
    typedef typename GridType::PositionType PositionType;
    typedef typename GridType::IndexType IndexType;
    typedef typename GridType::ValueType ValueType;
    typedef typename ScalarType<PositionType>::Type ScalarPositionType;

    FieldInterpolatorCICBase(const GridType& grid) :
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
    const GridType& grid;
    PositionType inverseStep;
    PositionType normalizedOrigin, normalizedOriginStaggered;

};

} // namespace pica::internal


template<typename Real>
class FieldInterpolatorCIC<YeeGrid<One, Real> > : public internal::FieldInterpolatorCICBase<One, Real> {
public:

    FieldInterpolatorCIC(const GridType& grid) : internal::FieldInterpolatorCICBase<One, Real>(grid)
    {
    }

    void get(const PositionType position, Vector3<ValueType>& e, Vector3<ValueType>& b)
    {
        IndexType indexCollocated, indexStaggered;
        PositionType coeffCollocated, coeffStaggered;
        getIndexCoeff(position, indexCollocated, coeffCollocated, indexStaggered, coeffStaggered);
        e.x = interpolate(&GridType::ex, IndexType(indexCollocated.x), PositionType(coeffCollocated.x));
        e.y = interpolate(&GridType::ey, IndexType(indexStaggered.x), PositionType(coeffStaggered.x));
        e.z = interpolate(&GridType::ey, IndexType(indexStaggered.x), PositionType(coeffStaggered.x));
        b.x = interpolate(&GridType::bx, IndexType(indexStaggered.x), PositionType(coeffStaggered.x));
        b.y = interpolate(&GridType::by, IndexType(indexCollocated.x), PositionType(coeffCollocated.x));
        b.z = interpolate(&GridType::bz, IndexType(indexCollocated.x), PositionType(coeffCollocated.x));
    }

private:

    typedef ValueType(GridType::*FieldComponent1d)(int) const;
    ValueType interpolate(FieldComponent1d component, IndexType baseIndex, PositionType coeff) const
    {
        return (grid.*component)(baseIndex.x) * (1.0 - coeff.x) +
            (grid.*component)(baseIndex.x + 1) * coeff.x;
    }

};


template<typename Real>
class FieldInterpolatorCIC<YeeGrid<Two, Real> > : public internal::FieldInterpolatorCICBase<Two, Real> {
public:

    FieldInterpolatorCIC(const GridType& grid) : internal::FieldInterpolatorCICBase<Two, Real>(grid)
    {
    }

    void get(const PositionType position, Vector3<ValueType>& e, Vector3<ValueType>& b)
    {
        IndexType indexCollocated, indexStaggered;
        PositionType coeffCollocated, coeffStaggered;
        getIndexCoeff(position, indexCollocated, coeffCollocated, indexStaggered, coeffStaggered);
        e.x = interpolate(&GridType::ex, IndexType(indexCollocated.x, indexStaggered.y),
            PositionType(coeffCollocated.x, coeffStaggered.y));
        e.y = interpolate(&GridType::ey, IndexType(indexStaggered.x, indexCollocated.y),
            PositionType(coeffStaggered.x, coeffCollocated.y));
        e.z = interpolate(&GridType::ey, IndexType(indexStaggered.x, indexStaggered.y),
            PositionType(coeffStaggered.x, coeffStaggered.y));
        b.x = interpolate(&GridType::bx, IndexType(indexStaggered.x, indexCollocated.y),
            PositionType(coeffStaggered.x, coeffCollocated.y));
        b.y = interpolate(&GridType::by, IndexType(indexCollocated.x, indexStaggered.y),
            PositionType(coeffCollocated.x, coeffStaggered.y));
        b.z = interpolate(&GridType::bz, IndexType(indexCollocated.x, indexCollocated.y),
            PositionType(coeffCollocated.x, coeffCollocated.y));
    }

private:

    typedef ValueType(GridType::*FieldComponent2d)(int, int) const;
    ValueType interpolate(FieldComponent2d component, const IndexType& baseIndex, const PositionType& coeff) const
    {
        return
            (grid.*component)(baseIndex.x, baseIndex.y) * (1.0 - coeff.x) * (1.0 - coeff.y) +
            (grid.*component)(baseIndex.x, baseIndex.y + 1) * (1.0 - coeff.x) * coeff.y +
            (grid.*component)(baseIndex.x + 1, baseIndex.y) * coeff.x * (1.0 - coeff.y) +
            (grid.*component)(baseIndex.x + 1, baseIndex.y + 1) * coeff.x * coeff.y;
    }

};


template<typename Real>
class FieldInterpolatorCIC<YeeGrid<Three, Real> > : public internal::FieldInterpolatorCICBase<Three, Real> {
public:

    FieldInterpolatorCIC(const GridType& grid) : internal::FieldInterpolatorCICBase<Three, Real>(grid)
    {
    }

    void get(const PositionType position, Vector3<ValueType>& e, Vector3<ValueType>& b)
    {
        IndexType indexCollocated, indexStaggered;
        PositionType coeffCollocated, coeffStaggered;
        getIndexCoeff(position, indexCollocated, coeffCollocated, indexStaggered, coeffStaggered);
        e.x = interpolate(&GridType::ex, IndexType(indexCollocated.x, indexStaggered.y, indexStaggered.z),
            PositionType(coeffCollocated.x, coeffStaggered.y, coeffStaggered.z));
        e.y = interpolate(&GridType::ey, IndexType(indexStaggered.x, indexCollocated.y, indexStaggered.z),
            PositionType(coeffStaggered.x, coeffCollocated.y, coeffStaggered.z));
        e.z = interpolate(&GridType::ey, IndexType(indexStaggered.x, indexStaggered.y, indexCollocated.z),
            PositionType(coeffStaggered.x, coeffStaggered.y, coeffCollocated.z));
        b.x = interpolate(&GridType::bx, IndexType(indexStaggered.x, indexCollocated.y, indexCollocated.z),
            PositionType(coeffStaggered.x, coeffCollocated.y, coeffCollocated.z));
        b.y = interpolate(&GridType::by, IndexType(indexCollocated.x, indexStaggered.y, indexCollocated.z),
            PositionType(coeffCollocated.x, coeffStaggered.y, coeffCollocated.z));
        b.z = interpolate(&GridType::bz, IndexType(indexCollocated.x, indexCollocated.y, indexStaggered.z),
            PositionType(coeffCollocated.x, coeffCollocated.y, coeffStaggered.z));
    }

private:
    
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

};


} // namespace pica


#endif
