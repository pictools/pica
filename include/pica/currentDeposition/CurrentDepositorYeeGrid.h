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
    typedef typename GridType::ArrayType ArrayType;
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

    using typename internal::CurrentDepositorYeeGridCICBase<One, Real>::GridType;
    using typename internal::CurrentDepositorYeeGridCICBase<One, Real>::IndexType;
    using typename internal::CurrentDepositorYeeGridCICBase<One, Real>::PositionType;
    using typename internal::CurrentDepositorYeeGridCICBase<One, Real>::ValueType;
    using typename internal::CurrentDepositorYeeGridCICBase<One, Real>::ArrayType;

    CurrentDepositorCIC(GridType& grid) : internal::CurrentDepositorYeeGridCICBase<One, Real>(grid)
    {}
 
    void deposit(const PositionType& position, const Vector3<ValueType>& current)
    {
        IndexType indexCollocated, indexStaggered;
        PositionType coeffCollocated, coeffStaggered;
        this->getIndexCoeff(position, indexCollocated, coeffCollocated, indexStaggered, coeffStaggered);
        deposit(this->grid.jx(), IndexType(indexCollocated.x), PositionType(coeffCollocated.x), current.x);
        deposit(this->grid.jy(), IndexType(indexStaggered.x), PositionType(coeffStaggered.x), current.y);
        deposit(this->grid.jz(), IndexType(indexStaggered.x), PositionType(coeffStaggered.x), current.z);
    }

private:

    void deposit(ArrayType& component, IndexType baseIndex, PositionType coeff, Real value) const
    {
        component(baseIndex.x) += (1.0 - coeff.x) * value;
        component(baseIndex.x + 1) += coeff.x * value;
    }

};


template<typename Real>
class CurrentDepositorCIC<YeeGrid<Two, Real> > : public internal::CurrentDepositorYeeGridCICBase<Two, Real> {
public:

    using typename internal::CurrentDepositorYeeGridCICBase<Two, Real>::GridType;
    using typename internal::CurrentDepositorYeeGridCICBase<Two, Real>::IndexType;
    using typename internal::CurrentDepositorYeeGridCICBase<Two, Real>::PositionType;
    using typename internal::CurrentDepositorYeeGridCICBase<Two, Real>::ValueType;
    using typename internal::CurrentDepositorYeeGridCICBase<Two, Real>::ArrayType;

    CurrentDepositorCIC(GridType& grid) : internal::CurrentDepositorYeeGridCICBase<Two, Real>(grid)
    {}

    void deposit(const PositionType& position, const Vector3<ValueType>& current)
    {
        IndexType indexCollocated, indexStaggered;
        PositionType coeffCollocated, coeffStaggered;
        this->getIndexCoeff(position, indexCollocated, coeffCollocated, indexStaggered, coeffStaggered);
        deposit(this->grid.jx(), IndexType(indexCollocated.x, indexStaggered.y),
            PositionType(coeffCollocated.x, coeffStaggered.y), current.x);
        deposit(this->grid.jy(), IndexType(indexStaggered.x, indexCollocated.y),
            PositionType(coeffStaggered.x, coeffCollocated.y), current.y);
        deposit(this->grid.jz(), IndexType(indexStaggered.x, indexStaggered.y),
            PositionType(coeffStaggered.x, coeffStaggered.y), current.z);
    }

private:

    void deposit(ArrayType component, IndexType baseIndex, PositionType coeff, Real value) const
    {
        component(baseIndex.x, baseIndex.y) += (1.0 - coeff.x) * (1.0 - coeff.y) * value;
        component(baseIndex.x, baseIndex.y + 1) += (1.0 - coeff.x) * coeff.y * value;
        component(baseIndex.x + 1, baseIndex.y) += coeff.x * (1.0 - coeff.y) * value;
        component(baseIndex.x + 1, baseIndex.y + 1) += coeff.x * coeff.y * value;
    }

};


template<typename Real>
class CurrentDepositorCIC<YeeGrid<Three, Real> > : public internal::CurrentDepositorYeeGridCICBase<Three, Real> {
public:

    using typename internal::CurrentDepositorYeeGridCICBase<Three, Real>::GridType;
    using typename internal::CurrentDepositorYeeGridCICBase<Three, Real>::IndexType;
    using typename internal::CurrentDepositorYeeGridCICBase<Three, Real>::PositionType;
    using typename internal::CurrentDepositorYeeGridCICBase<Three, Real>::ValueType;
    using typename internal::CurrentDepositorYeeGridCICBase<Three, Real>::ArrayType;

    CurrentDepositorCIC(GridType& grid) : internal::CurrentDepositorYeeGridCICBase<Three, Real>(grid)
    {}

    void deposit(const PositionType& position, const Vector3<ValueType>& current)
    {
        IndexType indexCollocated, indexStaggered;
        PositionType coeffCollocated, coeffStaggered;
        this->getIndexCoeff(position, indexCollocated, coeffCollocated, indexStaggered, coeffStaggered);
        deposit(this->grid.jx(), IndexType(indexCollocated.x, indexStaggered.y, indexStaggered.z),
            PositionType(coeffCollocated.x, coeffStaggered.y, coeffStaggered.z), current.x);
        deposit(this->grid.jy(), IndexType(indexStaggered.x, indexCollocated.y, indexStaggered.z),
            PositionType(coeffStaggered.x, coeffCollocated.y, coeffStaggered.z), current.y);
        deposit(this->grid.jz(), IndexType(indexStaggered.x, indexStaggered.y, indexCollocated.z),
            PositionType(coeffStaggered.x, coeffStaggered.y, coeffCollocated.z), current.z);
    }

private:

    void deposit(ArrayType& component, IndexType baseIndex, PositionType coeff, Real value) const
    {
        component(baseIndex.x, baseIndex.y, baseIndex.z) += (1.0 - coeff.x) * (1.0 - coeff.y) * (1.0 - coeff.z) * value;
        component(baseIndex.x, baseIndex.y, baseIndex.z + 1) += (1.0 - coeff.x) * (1.0 - coeff.y) * coeff.z * value;
        component(baseIndex.x, baseIndex.y + 1, baseIndex.z) += (1.0 - coeff.x) * coeff.y * (1.0 - coeff.z) * value;
        component(baseIndex.x, baseIndex.y + 1, baseIndex.z + 1) += (1.0 - coeff.x) * coeff.y * coeff.z * value;
        component(baseIndex.x + 1, baseIndex.y, baseIndex.z) += coeff.x * (1.0 - coeff.y) * (1.0 - coeff.z) * value;
        component(baseIndex.x + 1, baseIndex.y, baseIndex.z + 1) += coeff.x * (1.0 - coeff.y) * coeff.z * value;
        component(baseIndex.x + 1, baseIndex.y + 1, baseIndex.z) += coeff.x * coeff.y * (1.0 - coeff.z) * value;
        component(baseIndex.x + 1, baseIndex.y + 1, baseIndex.z + 1) += coeff.x * coeff.y * coeff.z * value;
    }

};


template<typename Real>
class CurrentDepositorCICSupercell {
public:

    typedef YeeGrid<Three, Real> GridType;
    typedef typename GridType::PositionType PositionType;
    typedef typename GridType::IndexType IndexType;
    typedef typename GridType::ValueType ValueType;
    typedef typename ScalarType<PositionType>::Type ScalarPositionType;

    CurrentDepositorCICSupercell(GridType& grid, const PositionType& minPosition,
        const IndexType& numCellsPerSupercell) :
        grid(grid),
        inverseStep(inverse(grid.getStep())),
        normalizedOrigin((minPosition - grid.getStep()) / grid.getStep()),
        normalizedOriginStaggered(normalizedOrigin + ones<Three, ScalarPositionType>() * static_cast<ScalarPositionType>(0.5))
    {
        baseIdx = truncate(minPosition / grid.getStep()) + Int3(1, 1, 1);
        size = numCellsPerSupercell + IndexType(3, 3, 3);
        /// Poor man's check
        for (int d = 0; d < 3; d++)
            if (size[d] > maxSize)
                std::cout << "ERROR: too large supercell size for depositor\n";

        for (int i = 0; i < size.x; i++)
        for (int j = 0; j < size.y; j++)
        for (int k = 0; k < size.z; k++) {
            jx[i][j][k] = 0.0;
            jy[i][j][k] = 0.0;
            jz[i][j][k] = 0.0;
        }
    }

    ~CurrentDepositorCICSupercell()
    {
        double normalization = 1.0 / grid.getStep().volume();
        for (int i = 0; i < size.x; i++)
        for (int j = 0; j < size.y; j++)
        for (int k = 0; k < size.z; k++)
        {
            IndexType gridxIdx = baseIdx + IndexType(i, j, k);
            grid.jx(gridxIdx.x, gridxIdx.y, gridxIdx.z) += jx[i][j][k] * normalization;
            grid.jy(gridxIdx.x, gridxIdx.y, gridxIdx.z) += jy[i][j][k] * normalization;
            grid.jz(gridxIdx.x, gridxIdx.y, gridxIdx.z) += jz[i][j][k] * normalization;
        }
    }

    void deposit(const PositionType& position, const Vector3<ValueType>& current)
    {
        IndexType indexCollocated, indexStaggered;
        PositionType coeffCollocated, coeffStaggered;
        this->getIndexCoeff(position, indexCollocated, coeffCollocated, indexStaggered, coeffStaggered);
        deposit(jx, IndexType(indexCollocated.x, indexStaggered.y, indexStaggered.z),
            PositionType(coeffCollocated.x, coeffStaggered.y, coeffStaggered.z), current.x);
        deposit(jy, IndexType(indexStaggered.x, indexCollocated.y, indexStaggered.z),
            PositionType(coeffStaggered.x, coeffCollocated.y, coeffStaggered.z), current.y);
        deposit(jz, IndexType(indexStaggered.x, indexStaggered.y, indexCollocated.z),
            PositionType(coeffStaggered.x, coeffStaggered.y, coeffCollocated.z), current.z);
    }

private:

    GridType& grid;
    PositionType inverseStep;
    PositionType normalizedOrigin, normalizedOriginStaggered;
    IndexType baseIdx, size;
    static const int maxSize = 7; // corresponds to max supercell size 4
    ValueType jx[maxSize][maxSize][maxSize];
    ValueType jy[maxSize][maxSize][maxSize];
    ValueType jz[maxSize][maxSize][maxSize];

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

    void deposit(ValueType component[maxSize][maxSize][maxSize], IndexType baseIndex,
        PositionType coeff, Real value) const
    {
        component[baseIndex.x][baseIndex.y][baseIndex.z] += (1.0 - coeff.x) * (1.0 - coeff.y) * (1.0 - coeff.z) * value;
        component[baseIndex.x][baseIndex.y][baseIndex.z + 1] += (1.0 - coeff.x) * (1.0 - coeff.y) * coeff.z * value;
        component[baseIndex.x][baseIndex.y + 1][baseIndex.z] += (1.0 - coeff.x) * coeff.y * (1.0 - coeff.z) * value;
        component[baseIndex.x][baseIndex.y + 1][baseIndex.z + 1] += (1.0 - coeff.x) * coeff.y * coeff.z * value;
        component[baseIndex.x + 1][baseIndex.y][baseIndex.z] += coeff.x * (1.0 - coeff.y) * (1.0 - coeff.z) * value;
        component[baseIndex.x + 1][baseIndex.y][baseIndex.z + 1] += coeff.x * (1.0 - coeff.y) * coeff.z * value;
        component[baseIndex.x + 1][baseIndex.y + 1][baseIndex.z] += coeff.x * coeff.y * (1.0 - coeff.z) * value;
        component[baseIndex.x + 1][baseIndex.y + 1][baseIndex.z + 1] += coeff.x * coeff.y * coeff.z * value;
    }

};


} // namespace pica


#endif
