#ifndef PICA_YEE_GRID_H
#define PICA_YEE_GRID_H


#include "pica/grid/Grid.h"


namespace pica {


// In 3D index(i, j, k) refers to cell(i, j, k) being the cell between nodes(i, j, k)
//        and (i + 1, j + 1, k + 1).It corresponds to Ex(i, j + 1 / 2, k + 1 / 2),
//        Ey(i + 1 / 2, j, k + 1 / 2), Ez(i + 1 / 2, j + 1 / 2, k), same for the components of J,
//        Bx(i + 1 / 2, j, k), By(i, j + 1 / 2, k), Bz(i, j, k + 1 / 2).
template<Dimension dimension, typename Real = double>
class YeeGrid: public Grid_<dimension, Real> {
public:
    YeeGrid(const IndexType& size) :
        Grid_(size) {}
};


} // namespace pica


#endif
