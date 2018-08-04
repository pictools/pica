#ifndef PICA_FIELD_INTERPOLATOR_YEE_GRID_H
#define PICA_FIELD_INTERPOLATOR_YEE_GRID_H


#include "pica/grid/YeeGrid.h"


namespace pica {


namespace internal {

template<Dimension dimension, typename Real>
class FieldInterpolatorYeeGridCICBase {
public:
    typedef YeeGrid<dimension, Real> GridType;
    typedef typename GridType::PositionType PositionType;
    typedef typename GridType::IndexType IndexType;
    typedef typename GridType::ValueType ValueType;
    typedef typename GridType::ArrayType ArrayType;
    typedef typename ScalarType<PositionType>::Type ScalarPositionType;

    FieldInterpolatorYeeGridCICBase(const GridType& grid) :
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
class FieldInterpolatorCIC<YeeGrid<One, Real> > : public internal::FieldInterpolatorYeeGridCICBase<One, Real> {
public:

    using typename internal::FieldInterpolatorYeeGridCICBase<One, Real>::GridType;
    using typename internal::FieldInterpolatorYeeGridCICBase<One, Real>::IndexType;
    using typename internal::FieldInterpolatorYeeGridCICBase<One, Real>::PositionType;
    using typename internal::FieldInterpolatorYeeGridCICBase<One, Real>::ValueType;
    using typename internal::FieldInterpolatorYeeGridCICBase<One, Real>::ArrayType;

    FieldInterpolatorCIC(const GridType& grid) : internal::FieldInterpolatorYeeGridCICBase<One, Real>(grid)
    {
    }

    void get(const PositionType position, Vector3<ValueType>& e, Vector3<ValueType>& b)
    {
        IndexType indexCollocated, indexStaggered;
        PositionType coeffCollocated, coeffStaggered;
        this->getIndexCoeff(position, indexCollocated, coeffCollocated, indexStaggered, coeffStaggered);
        e.x = interpolate(this->grid.ex(), IndexType(indexCollocated.x), PositionType(coeffCollocated.x));
        e.y = interpolate(this->grid.ey(), IndexType(indexStaggered.x), PositionType(coeffStaggered.x));
        e.z = interpolate(this->grid.ez(), IndexType(indexStaggered.x), PositionType(coeffStaggered.x));
        b.x = interpolate(this->grid.bx(), IndexType(indexStaggered.x), PositionType(coeffStaggered.x));
        b.y = interpolate(this->grid.by(), IndexType(indexCollocated.x), PositionType(coeffCollocated.x));
        b.z = interpolate(this->grid.bz(), IndexType(indexCollocated.x), PositionType(coeffCollocated.x));
    }

private:

    ValueType interpolate(const ArrayType& component, IndexType baseIndex, PositionType coeff) const
    {
        return component(baseIndex.x) * (1.0 - coeff.x) + component(baseIndex.x + 1) * coeff.x;
    }

};


template<typename Real>
class FieldInterpolatorCIC<YeeGrid<Two, Real> > : public internal::FieldInterpolatorYeeGridCICBase<Two, Real> {
public:

    using typename internal::FieldInterpolatorYeeGridCICBase<Two, Real>::GridType;
    using typename internal::FieldInterpolatorYeeGridCICBase<Two, Real>::IndexType;
    using typename internal::FieldInterpolatorYeeGridCICBase<Two, Real>::PositionType;
    using typename internal::FieldInterpolatorYeeGridCICBase<Two, Real>::ValueType;
    using typename internal::FieldInterpolatorYeeGridCICBase<Two, Real>::ArrayType;

    FieldInterpolatorCIC(const GridType& grid) : internal::FieldInterpolatorYeeGridCICBase<Two, Real>(grid)
    {
    }

    void get(const PositionType position, Vector3<ValueType>& e, Vector3<ValueType>& b)
    {
        IndexType indexCollocated, indexStaggered;
        PositionType coeffCollocated, coeffStaggered;
        this->getIndexCoeff(position, indexCollocated, coeffCollocated, indexStaggered, coeffStaggered);
        e.x = interpolate(this->grid.ex(), IndexType(indexCollocated.x, indexStaggered.y),
            PositionType(coeffCollocated.x, coeffStaggered.y, coeffStaggered.z));
        e.y = interpolate(this->grid.ey(), IndexType(indexStaggered.x, indexCollocated.y),
            PositionType(coeffStaggered.x, coeffCollocated.y, coeffStaggered.z));
        e.z = interpolate(this->grid.ez(), IndexType(indexStaggered.x, indexStaggered.y),
            PositionType(coeffStaggered.x, coeffStaggered.y, coeffCollocated.z));
        b.x = interpolate(this->grid.bx(), IndexType(indexStaggered.x, indexCollocated.y),
            PositionType(coeffStaggered.x, coeffCollocated.y, coeffCollocated.z));
        b.y = interpolate(this->grid.by(), IndexType(indexCollocated.x, indexStaggered.y),
            PositionType(coeffCollocated.x, coeffStaggered.y, coeffCollocated.z));
        b.z = interpolate(this->grid.bz(), IndexType(indexCollocated.x, indexCollocated.y),
            PositionType(coeffCollocated.x, coeffCollocated.y, coeffStaggered.z));

    }

private:

    ValueType interpolate(const ArrayType component, const IndexType& baseIndex, const PositionType& coeff) const
    {
        return
            component(baseIndex.x, baseIndex.y) * (1.0 - coeff.x) * (1.0 - coeff.y) +
            component(baseIndex.x, baseIndex.y + 1) * (1.0 - coeff.x) * coeff.y +
            component(baseIndex.x + 1, baseIndex.y) * coeff.x * (1.0 - coeff.y) +
            component(baseIndex.x + 1, baseIndex.y + 1) * coeff.x * coeff.y;
    }

};


template<typename Real>
class FieldInterpolatorCIC<YeeGrid<Three, Real> > : public internal::FieldInterpolatorYeeGridCICBase<Three, Real> {
public:

    using typename internal::FieldInterpolatorYeeGridCICBase<Three, Real>::GridType;
    using typename internal::FieldInterpolatorYeeGridCICBase<Three, Real>::IndexType;
    using typename internal::FieldInterpolatorYeeGridCICBase<Three, Real>::PositionType;
    using typename internal::FieldInterpolatorYeeGridCICBase<Three, Real>::ValueType;
    using typename internal::FieldInterpolatorYeeGridCICBase<Three, Real>::ArrayType;

    FieldInterpolatorCIC(const GridType& grid) : internal::FieldInterpolatorYeeGridCICBase<Three, Real>(grid)
    {
    }

    void get(const PositionType position, Vector3<ValueType>& e, Vector3<ValueType>& b)
    {
        IndexType indexCollocated, indexStaggered;
        PositionType coeffCollocated, coeffStaggered;
        this->getIndexCoeff(position, indexCollocated, coeffCollocated, indexStaggered, coeffStaggered);
        e.x = interpolate(this->grid.ex(), IndexType(indexCollocated.x, indexStaggered.y, indexStaggered.z),
            PositionType(coeffCollocated.x, coeffStaggered.y, coeffStaggered.z));
        e.y = interpolate(this->grid.ey(), IndexType(indexStaggered.x, indexCollocated.y, indexStaggered.z),
            PositionType(coeffStaggered.x, coeffCollocated.y, coeffStaggered.z));
        e.z = interpolate(this->grid.ez(), IndexType(indexStaggered.x, indexStaggered.y, indexCollocated.z),
            PositionType(coeffStaggered.x, coeffStaggered.y, coeffCollocated.z));
        b.x = interpolate(this->grid.bx(), IndexType(indexStaggered.x, indexCollocated.y, indexCollocated.z),
            PositionType(coeffStaggered.x, coeffCollocated.y, coeffCollocated.z));
        b.y = interpolate(this->grid.by(), IndexType(indexCollocated.x, indexStaggered.y, indexCollocated.z),
            PositionType(coeffCollocated.x, coeffStaggered.y, coeffCollocated.z));
        b.z = interpolate(this->grid.bz(), IndexType(indexCollocated.x, indexCollocated.y, indexStaggered.z),
            PositionType(coeffCollocated.x, coeffCollocated.y, coeffStaggered.z));
    }

private:

    ValueType interpolate(const ArrayType& component, const IndexType& baseIndex, const PositionType& coeff) const
    {
        return
            component(baseIndex.x, baseIndex.y, baseIndex.z) * (1.0 - coeff.x) * (1.0 - coeff.y) * (1.0 - coeff.z) +
            component(baseIndex.x, baseIndex.y, baseIndex.z + 1) * (1.0 - coeff.x) * (1.0 - coeff.y) * coeff.z +
            component(baseIndex.x, baseIndex.y + 1, baseIndex.z) * (1.0 - coeff.x) * coeff.y * (1.0 - coeff.z) +
            component(baseIndex.x, baseIndex.y + 1, baseIndex.z + 1) * (1.0 - coeff.x) * coeff.y * coeff.z +
            component(baseIndex.x + 1, baseIndex.y, baseIndex.z) * coeff.x * (1.0 - coeff.y) * (1.0 - coeff.z) +
            component(baseIndex.x + 1, baseIndex.y, baseIndex.z + 1) * coeff.x * (1.0 - coeff.y) * coeff.z +
            component(baseIndex.x + 1, baseIndex.y + 1, baseIndex.z) * coeff.x * coeff.y * (1.0 - coeff.z) +
            component(baseIndex.x + 1, baseIndex.y + 1, baseIndex.z + 1) * coeff.x * coeff.y * coeff.z;
    }

};


template<typename Real>
class FieldInterpolatorCICSupercell {
public:

    typedef YeeGrid<Three, Real> GridType;
    typedef typename GridType::PositionType PositionType;
    typedef typename GridType::IndexType IndexType;
    typedef typename GridType::ValueType ValueType;
    typedef typename ScalarType<PositionType>::Type ScalarPositionType;

    FieldInterpolatorCICSupercell(const GridType& grid, const PositionType& minPosition,
        const IndexType& numCellsPerSupercell) :
        grid(grid),
        inverseStep(inverse(grid.getStep())),
        normalizedOrigin((minPosition - grid.getStep()) / grid.getStep()),
        normalizedOriginStaggered(normalizedOrigin + ones<Three, ScalarPositionType>() * static_cast<ScalarPositionType>(0.5))
    {
        IndexType size = numCellsPerSupercell + IndexType(2, 2, 2);
        IndexType baseIdx = truncate(minPosition / grid.getStep()) + Int3(1, 1, 1);
        /// Poor man's check
        for (int d = 0; d < 3; d++)
            if (size[d] > maxSize)
                std::cout << "ERROR: too large supercell size for interpolator\n";

        for (int i = 0; i < size.x; i++)
        for (int j = 0; j < size.y; j++)
        for (int k = 0; k < size.z; k++)
        {
            IndexType gridxIdx = baseIdx + IndexType(i, j, k);
            ex[i][j][k] = grid.ex(gridxIdx.x, gridxIdx.y, gridxIdx.z);
            ey[i][j][k] = grid.ey(gridxIdx.x, gridxIdx.y, gridxIdx.z);
            ez[i][j][k] = grid.ez(gridxIdx.x, gridxIdx.y, gridxIdx.z);
            bx[i][j][k] = grid.bx(gridxIdx.x, gridxIdx.y, gridxIdx.z);
            by[i][j][k] = grid.by(gridxIdx.x, gridxIdx.y, gridxIdx.z);
            bz[i][j][k] = grid.bz(gridxIdx.x, gridxIdx.y, gridxIdx.z); 
        }
    }

    void get(const PositionType position, Vector3<ValueType>& e, Vector3<ValueType>& b)
    {
        IndexType indexCollocated, indexStaggered;
        PositionType coeffCollocated, coeffStaggered;
        this->getIndexCoeff(position, indexCollocated, coeffCollocated, indexStaggered, coeffStaggered);
        e.x = interpolate(ex, IndexType(indexCollocated.x, indexStaggered.y, indexStaggered.z),
            PositionType(coeffCollocated.x, coeffStaggered.y, coeffStaggered.z));
        e.y = interpolate(ey, IndexType(indexStaggered.x, indexCollocated.y, indexStaggered.z),
            PositionType(coeffStaggered.x, coeffCollocated.y, coeffStaggered.z));
        e.z = interpolate(ez, IndexType(indexStaggered.x, indexStaggered.y, indexCollocated.z),
            PositionType(coeffStaggered.x, coeffStaggered.y, coeffCollocated.z));
        b.x = interpolate(bx, IndexType(indexStaggered.x, indexCollocated.y, indexCollocated.z),
            PositionType(coeffStaggered.x, coeffCollocated.y, coeffCollocated.z));
        b.y = interpolate(by, IndexType(indexCollocated.x, indexStaggered.y, indexCollocated.z),
            PositionType(coeffCollocated.x, coeffStaggered.y, coeffCollocated.z));
        b.z = interpolate(bz, IndexType(indexCollocated.x, indexCollocated.y, indexStaggered.z),
            PositionType(coeffCollocated.x, coeffCollocated.y, coeffStaggered.z));
    }

private:

    const GridType& grid;
    PositionType inverseStep;
    PositionType normalizedOrigin, normalizedOriginStaggered;
    static const int maxSize = 6; // corresponds to max supercell size 4
    ValueType ex[maxSize][maxSize][maxSize];
    ValueType ey[maxSize][maxSize][maxSize];
    ValueType ez[maxSize][maxSize][maxSize];
    ValueType bx[maxSize][maxSize][maxSize];
    ValueType by[maxSize][maxSize][maxSize];
    ValueType bz[maxSize][maxSize][maxSize];

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

    ValueType interpolate(ValueType component[maxSize][maxSize][maxSize],
        const IndexType& baseIndex, const PositionType& coeff) const
    {
        return
            component[baseIndex.x][baseIndex.y][baseIndex.z] * (1.0 - coeff.x) * (1.0 - coeff.y) * (1.0 - coeff.z) +
            component[baseIndex.x][baseIndex.y][baseIndex.z + 1] * (1.0 - coeff.x) * (1.0 - coeff.y) * coeff.z +
            component[baseIndex.x][baseIndex.y + 1][baseIndex.z] * (1.0 - coeff.x) * coeff.y * (1.0 - coeff.z) +
            component[baseIndex.x][baseIndex.y + 1][baseIndex.z + 1] * (1.0 - coeff.x) * coeff.y * coeff.z +
            component[baseIndex.x + 1][baseIndex.y][baseIndex.z] * coeff.x * (1.0 - coeff.y) * (1.0 - coeff.z) +
            component[baseIndex.x + 1][baseIndex.y][baseIndex.z + 1] * coeff.x * (1.0 - coeff.y) * coeff.z +
            component[baseIndex.x + 1][baseIndex.y + 1][baseIndex.z] * coeff.x * coeff.y * (1.0 - coeff.z) +
            component[baseIndex.x + 1][baseIndex.y + 1][baseIndex.z + 1] * coeff.x * coeff.y * coeff.z;
    }

};



} // namespace pica


#endif
