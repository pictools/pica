#ifndef PICA_CURRENT_DEPOSITOR_INTERFACE_H
#define PICA_CURRENT_DEPOSITOR_INTERFACE_H


#include "pica/math/Vectors.h"


namespace pica {


template<class Grid>
class CurrentDepositorCIC {
public:
    CurrentDepositorCIC(Grid& grid);
    void deposit(const typename Grid::PositionType& position, const Vector3<typename Grid::ValueType>& current);
};


} // namespace pica


#endif
