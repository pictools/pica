#ifndef PICA_FIELD_INTERPOLATOR_INTERFACE_H
#define PICA_FIELD_INTERPOLATOR_INTERFACE_H


#include "pica/math/Vectors.h"


namespace pica {


template<class Grid>
class FieldInterpolatorCIC {
public:
    FieldInterpolatorCIC(const Grid& grid);
    void get(const typename Grid::PositionType& position, Vector3<typename Grid::ValueType>& e, Vector3<typename Grid::ValueType>& b);
};


} // namespace pica


#endif
