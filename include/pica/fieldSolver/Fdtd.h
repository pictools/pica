#ifndef PICA_FDTD_H
#define PICA_FDTD_H


#include "pica/fieldSolver/FieldSolver.h"
#include "pica/grid/Grid.h"


namespace pica {

class Fdtd : public FieldSolver
{

public:

    Fdtd(const Parameters& parameters, Grid& grid);

    virtual void updateHalfB();
    virtual void updateE();

private:
    
    void updateHalfB3D();
    void updateHalfB2D();
    void updateHalfB1D();
    void updateE3D();
    void updateE2D();
    void updateE1D();

};

} // namespace pica

#endif
