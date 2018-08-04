#ifndef PICA_YEESOLVER_H
#define PICA_YEESOLVER_H


#include "pica/grid/YeeGrid.h"
#include "pica/threading/OpenMPHelper.h"


namespace pica {

class YeeSolver {
public:

    template<Dimension dimension, typename Real>
    void updateE(YeeGrid<dimension, Real>& grid, Real dt);

    template<Dimension dimension, typename Real>
    void updateB(YeeGrid<dimension, Real>& grid, Real dt);

private:

    template<Dimension dimension, typename Real>
    struct Implementation {};

};


template<Dimension dimension, typename Real>
void YeeSolver::updateE(YeeGrid<dimension, Real>& grid, Real dt)
{
    Implementation<dimension, Real>::updateE(grid, dt);
}

template<Dimension dimension, typename Real>
void YeeSolver::updateB(YeeGrid<dimension, Real>& grid, Real dt)
{
    Implementation<dimension, Real>::updateB(grid, dt);
}


template<typename Real>
struct YeeSolver::Implementation<One, Real> {
    static void updateE(YeeGrid<One, Real>& grid, Real dt)
    {
        typedef typename YeeGrid<One, Real>::ValueType ValueType;
        const ValueType coeffCurrent = -static_cast<ValueType>(4) * Constants<ValueType>::pi() * dt;
        typedef typename YeeGrid<One, Real>::PositionType PositionType;
        const ValueType cdt = Constants<ValueType>::c() * dt;
        const PositionType coeff = PositionType(cdt) / grid.getStep();
        typedef typename YeeGrid<One, Real>::IndexType IndexType;
        const IndexType begin = 0;
        const IndexType end = grid.getSize() - IndexType(1);
        #pragma omp parallel for
        for (int i = begin.x; i < end.x; i++) {
            grid.ex(i) += coeffCurrent * grid.jx(i);
            grid.ey(i) += coeffCurrent * grid.jy(i) - coeff.x * (grid.bz(i + 1) - grid.bz(i));
            grid.ez(i) += coeffCurrent * grid.jz(i) + coeff.x * (grid.by(i + 1) - grid.by(i));
        }
    }

    static void updateB(YeeGrid<One, Real>& grid, Real dt)
    {
        typedef typename YeeGrid<One, Real>::ValueType ValueType;
        typedef typename YeeGrid<One, Real>::PositionType PositionType;
        const ValueType cdt = Constants<ValueType>::c() * dt;
        const PositionType coeff = PositionType(cdt) / grid.getStep();
        typedef typename YeeGrid<One, Real>::IndexType IndexType;
        const IndexType begin = IndexType(1);
        const IndexType end = grid.getSize();
        #pragma omp parallel for
        for (int i = begin.x; i < end.x; i++) {
            grid.by(i) += coeff.x * (grid.ez(i) - grid.ez(i - 1));
            grid.bz(i) += -coeff.x * (grid.ey(i) - grid.ey(i - 1));
        }
    }
};

template<typename Real>
struct YeeSolver::Implementation<Two, Real> {
    static void updateE(YeeGrid<Two, Real>& grid, Real dt)
    {
        typedef typename YeeGrid<Two, Real>::ValueType ValueType;
        const ValueType coeffCurrent = -static_cast<ValueType>(4) * Constants<ValueType>::pi() * dt;
        typedef typename YeeGrid<Two, Real>::PositionType PositionType;
        const ValueType cdt = Constants<ValueType>::c() * dt;
        const PositionType coeff = PositionType(cdt, cdt) / grid.getStep();
        typedef typename YeeGrid<Two, Real>::IndexType IndexType;
        const IndexType begin(0, 0);
        const IndexType end = grid.getSize() - IndexType(1, 1);
        #pragma omp parallel for
        for (int i = begin.x; i < end.x; i++)
        for (int j = begin.y; j < end.y; j++) {
            grid.ex(i, j) += coeffCurrent * grid.jx(i, j) + coeff.y * (grid.bz(i, j + 1) - grid.bz(i, j));
            grid.ey(i, j) += coeffCurrent * grid.jy(i, j) - coeff.x * (grid.bz(i + 1, j) - grid.bz(i, j));
            grid.ez(i, j) += coeffCurrent * grid.jz(i, j) + coeff.x * (grid.by(i + 1, j) - grid.by(i, j)) -
                coeff.y * (grid.bx(i, j + 1) - grid.bx(i, j));
        }
        // Edges
        for (int i = begin.x; i < end.x; i++)
            grid.ey(i, end.y) += coeffCurrent * grid.jy(i, end.y) - coeff.x * (grid.bz(i + 1, end.y) - grid.bz(i, end.y));
        for (int j = begin.y; j < end.y; j++)
            grid.ex(end.x, j) += coeffCurrent * grid.jx(end.x, j) + coeff.y * (grid.bz(end.x, j + 1) - grid.bz(end.x, j));
    }

    static void updateB(YeeGrid<Two, Real>& grid, Real dt)
    {
        typedef typename YeeGrid<Two, Real>::ValueType ValueType;
        typedef typename YeeGrid<Two, Real>::PositionType PositionType;
        const ValueType cdt = Constants<ValueType>::c() * dt;
        const PositionType coeff = PositionType(cdt, cdt) / grid.getStep();
        typedef typename YeeGrid<Two, Real>::IndexType IndexType;
        const IndexType begin(1, 1);
        const IndexType end = grid.getSize();
        #pragma omp parallel for
        for (int i = begin.x; i < end.x; i++)
        for (int j = begin.y; j < end.y; j++) {
            grid.bx(i, j) += -coeff.y * (grid.ez(i, j) - grid.ez(i, j - 1));
            grid.by(i, j) += coeff.x * (grid.ez(i, j) - grid.ez(i - 1, j));
            grid.bz(i, j) += coeff.y * (grid.ex(i, j) - grid.ex(i, j - 1)) -
                coeff.x * (grid.ey(i, j) - grid.ey(i - 1, j));
        }
        // Edges
        for (int i = begin.x; i < end.x; i++)
            grid.by(i, 0) += coeff.x * (grid.ez(i, 0) - grid.ez(i - 1, 0));
        for (int j = begin.y; j < end.y; j++)
            grid.bx(0, j) += -coeff.y * (grid.ez(0, j) - grid.ez(0, j - 1));
    }
};

template<typename Real>
struct YeeSolver::Implementation<Three, Real> {
    static void updateE(YeeGrid<Three, Real>& grid, Real dt)
    {
        typedef typename YeeGrid<Three, Real>::ValueType ValueType;
        const ValueType coeffCurrent = -static_cast<ValueType>(4) * Constants<ValueType>::pi() * dt;
        typedef typename YeeGrid<Three, Real>::PositionType PositionType;
        const ValueType cdt = Constants<ValueType>::c() * dt;
        const PositionType coeff = PositionType(cdt, cdt, cdt) / grid.getStep();
        typedef typename YeeGrid<Three, Real>::IndexType IndexType;
        const IndexType begin(0, 0, 0);
        const IndexType end = grid.getSize() - IndexType(1, 1, 1);
        #pragma omp parallel for collapse(2)
        for (int i = begin.x; i < end.x; i++)
        for (int j = begin.y; j < end.y; j++)
        for (int k = begin.z; k < end.z; k++) {
            grid.ex(i, j, k) += coeffCurrent * grid.jx(i, j, k) +
                coeff.y * (grid.bz(i, j + 1, k) - grid.bz(i, j, k)) -
                coeff.z * (grid.by(i, j, k + 1) - grid.by(i, j, k));
            grid.ey(i, j, k) += coeffCurrent * grid.jy(i, j, k) +
                coeff.z * (grid.bx(i, j, k + 1) - grid.bx(i, j, k)) -
                coeff.x * (grid.bz(i + 1, j, k) - grid.bz(i, j, k));
            grid.ez(i, j, k) += coeffCurrent * grid.jz(i, j, k) +
                coeff.x * (grid.by(i + 1, j, k) - grid.by(i, j, k)) -
                coeff.y * (grid.bx(i, j + 1, k) - grid.bx(i, j, k));
        }
        // Edges
        #pragma omp parallel for
        for (int i = begin.x; i < end.x; i++)
        for (int j = begin.y; j < end.y; j++)
            grid.ez(i, j, end.z) += coeffCurrent * grid.jz(i, j, end.z) +
                coeff.x * (grid.by(i + 1, j, end.z) - grid.by(i, j, end.z)) -
                coeff.y * (grid.bx(i, j + 1, end.z) - grid.bx(i, j, end.z));
        #pragma omp parallel for
        for (int i = begin.x; i < end.x; i++)
        for (int k = begin.z; k < end.z; k++)
            grid.ey(i, end.y, k) += coeffCurrent * grid.jy(i, end.y, k) +
                coeff.z * (grid.bx(i, end.y, k + 1) - grid.bx(i, end.y, k)) -
                coeff.x * (grid.bz(i + 1, end.y, k) - grid.bz(i, end.y, k));
        #pragma omp parallel for
        for (int j = begin.y; j < end.y; j++)
        for (int k = begin.z; k < end.z; k++)
            grid.ex(end.x, j, k) += coeffCurrent * grid.jx(end.x, j, k) +
                coeff.y * (grid.bz(end.x, j + 1, k) - grid.bz(end.x, j, k)) -
                coeff.z * (grid.by(end.x, j, k + 1) - grid.by(end.x, j, k));
    }

    static void updateB(YeeGrid<Three, Real>& grid, Real dt)
    {
        typedef typename YeeGrid<Three, Real>::ValueType ValueType;
        typedef typename YeeGrid<Three, Real>::PositionType PositionType;
        const ValueType cdt = Constants<ValueType>::c() * dt;
        const PositionType coeff = PositionType(cdt, cdt, cdt) / grid.getStep();
        typedef typename YeeGrid<Three, Real>::IndexType IndexType;
        const IndexType begin(1, 1, 1);
        const IndexType end = grid.getSize();
        #pragma omp parallel for collapse(2)
        for (int i = begin.x; i < end.x; i++)
        for (int j = begin.y; j < end.y; j++)
        for (int k = begin.z; k < end.z; k++) {
            grid.bx(i, j, k) += coeff.z * (grid.ey(i, j, k) - grid.ey(i, j, k - 1)) -
                coeff.y * (grid.ez(i, j, k) - grid.ez(i, j - 1, k));
            grid.by(i, j, k) += coeff.x * (grid.ez(i, j, k) - grid.ez(i - 1, j, k)) -
                coeff.z * (grid.ex(i, j, k) - grid.ex(i, j, k - 1));
            grid.bz(i, j, k) += coeff.y * (grid.ex(i, j, k) - grid.ex(i, j - 1, k)) -
                coeff.x * (grid.ey(i, j, k) - grid.ey(i - 1, j, k));
        }
        // Edges
        #pragma omp parallel for
        for (int i = begin.x; i < end.x; i++)
        for (int j = begin.y; j < end.y; j++)
            grid.bz(i, j, 0) += coeff.y * (grid.ex(i, j, 0) - grid.ex(i, j - 1, 0)) -
                coeff.x * (grid.ey(i, j, 0) - grid.ey(i - 1, j, 0));
        #pragma omp parallel for
        for (int i = begin.x; i < end.x; i++)
        for (int k = begin.z; k < end.z; k++)
            grid.by(i, 0, k) += coeff.x * (grid.ez(i, 0, k) - grid.ez(i - 1, 0, k)) -
                coeff.z * (grid.ex(i, 0, k) - grid.ex(i, 0, k - 1));
        #pragma omp parallel for
        for (int j = begin.y; j < end.y; j++)
        for (int k = begin.z; k < end.z; k++)
            grid.bx(0, j, k) += coeff.z * (grid.ey(0, j, k) - grid.ey(0, j, k - 1)) -
                coeff.y * (grid.ez(0, j, k) - grid.ez(0, j - 1, k));
    }
};

} // namespace pica

#endif
