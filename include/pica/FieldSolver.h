#ifndef PICA_FIELDSOLVER_H
#define PICA_FIELDSOLVER_H


#include "pica/grid/Grid.h"
#include "Parameters.h"

#include <memory>


namespace pica {

// Base class for field solvers on Grid.
// The main method doStep uses template method pattern.
class FieldSolver
{

public:

    FieldSolver(const Parameters& parameters, Grid& grid);
 
    virtual void updateHalfB() = 0;
    virtual void updateE() = 0;

public:

    Parameters parameters;
    Grid* grid;
    // Index space being updated in form [begin, end).
    Int3 updateBAreaBegin, updateBAreaEnd;
    Int3 updateEAreaBegin, updateEAreaEnd;
    // Index space being updated in non-PML area.
    Int3 internalBAreaBegin, internalBAreaEnd;
    Int3 internalEAreaBegin, internalEAreaEnd;

    void updateDims();
    void updateInternalDims();
    void updateA();

private:

    // Copy and assignment are disallowed.
    FieldSolver(const FieldSolver &);
    FieldSolver & operator =(const FieldSolver &);
};

} // namespace pica

#endif

