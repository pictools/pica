#ifndef PICA_CONSTANTS_H
#define PICA_CONSTANTS_H


#include "pica/math/FP.h"


namespace pica {

/* Mathematical and physical constants in CGS. */
namespace constants
{
    const FP c = 29979245800.0;
    const FP LightVelocity = c;
    const FP electronCharge = -4.80320427e-10;
    const FP ElectronCharge = electronCharge;
    const FP electronMass = 9.10938215e-28;
    const FP ElectronMass = electronMass;
    const FP protonMass = 1.672622964e-24;
    const FP pi = 3.14159265358;
    const FP Pi = pi;
    const FP reducedPlanck = 1.0545716818e-27;
    const FP planckConstant = reducedPlanck;
    const FP electronVolt = 1.60217656535e-12;
    const FP eV = electronVolt;
    const FP MeV = 1e6 * eV;
} // namespace pica::constants

} // namespace pica


#endif
