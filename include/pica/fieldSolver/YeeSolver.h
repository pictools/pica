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
        const ValueType coeff = Constants<ValueType>::c() * dt / grid.getSteps();
        typedef typename YeeGrid<One, Real>::IndexType IndexType;
        const IndexType begin = 0;
        const IndexType end = grid.getSize() - 1;
        #pragma omp parallel for
        for (int i = begin; i < end; i++) {
            grid.ex(i) += coeffCurrent * grid.jx(i);
            grid.ey(i) += coeffCurrent * grid.jy(i) - coeff * (grid.bz(i + 1) - grid.bz(i));
            grid.ez(i) += coeffCurrent * grid.jz(i) + coeff * (grid.by(i + 1) - grid.by(i));
        }
    }

    static void updateB(YeeGrid<One, Real>& grid, Real dt)
    {
        typedef typename YeeGrid<One, Real>::ValueType ValueType;
        const ValueType coeff = Constants<ValueType>::c() * dt / grid.getSteps();
        typedef typename YeeGrid<One, Real>::IndexType IndexType;
        const IndexType begin = 1;
        const IndexType end = grid.getSize();
        #pragma omp parallel for
        for (int i = begin; i < end; i++) {
            grid.by(i) += coeff * (grid.ez(i) - grid.ez(i - 1));
            grid.bz(i) += -coeff * (grid.ey(i) - grid.ey(i - 1));
        }
    }
};

template<typename Real>
struct YeeSolver::Implementation<Two, Real> {
    static void updateE(YeeGrid<Two, Real>& grid, Real dt)
    {
        typedef typename YeeGrid<Two, Real>::ValueType ValueType;
        const ValueType coeffCurrent = -static_cast<ValueType>(4) * Constants<ValueType>::pi() * dt;
        typedef typename YeeGrid<Two, Real>::StepsType StepsType;
        const ValueType cdt = Constants<ValueType>::c() * dt;
        const StepsType coeff = StepsType(cdt, cdt) / grid.getSteps();
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
        typedef typename YeeGrid<Two, Real>::StepsType StepsType;
        const ValueType cdt = Constants<ValueType>::c() * dt;
        const StepsType coeff = StepsType(cdt, cdt) / grid.getSteps();
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
    }
};

template<typename Real>
struct YeeSolver::Implementation<Three, Real> {
    //static void updateE(YeeGrid<Three, Real>& grid, Real dt);
    //static void updateB(YeeGrid<Three, Real>& grid, Real dt);
};

} // namespace pica

#endif
