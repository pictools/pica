#ifndef PICA_PARTICLE_H
#define PICA_PARTICLE_H


#include "pica/math/Constants.h"
#include "pica/math/FP.h"
#include "pica/math/Vectors.h"

#include <cmath>
#include <stddef.h>
#include <string>
#include <vector>


namespace pica {


template<Dimension dimension>
class Particle {
public:

    // Types for conforming ParticleInterface
    typedef FP Real;
    typedef typename VectorTypeHelper<dimension, Real>::Type PositionType;
    typedef typename VectorTypeHelper<Three, Real>::Type MomentumType;
    typedef Real GammaType;
    typedef Real MassType;
    typedef Real ChargeType;
    typedef float FactorType;

    Particle() :
        factor(1),
        mass(0.0),
        charge(0.0),
        invGamma(1.0) {}

    Particle(const PositionType& position, const MomentumType& momentum,
        MassType mass, ChargeType charge, FactorType factor = 1) :
        position(position), mass(mass), charge(charge), factor(factor)
    {
        setMomentum(momentum);
    }

    PositionType getPosition() const { return position; }
    void setPosition(const PositionType& newPosition) { position = newPosition; }

    MomentumType getMomentum() const
    {
        return p * Constants<MassType>::c() * mass;
    }

    void setMomentum(const MomentumType& newMomentum)
    {
        p = newMomentum / (Constants<GammaType>::c() * mass);
        invGamma = static_cast<GammaType>(1.0) / sqrt(static_cast<GammaType>(1.0) + p.norm2());
    }

    void setP(const  MomentumType& newP)
    {
        p = newP;
        invGamma = static_cast<GammaType>(1.0) / sqrt(static_cast<GammaType>(1.0) + p.norm2());
    }

    MomentumType getP() const { return p; }

    MomentumType getVelocity() const
    {
        return p * (Constants<GammaType>::c() * invGamma);
    }

    void setVelocity(const MomentumType& newVelocity)
    {
        p = newVelocity / sqrt(constants::c * constants::c - newVelocity.norm2());
        invGamma = (FP)1 / sqrt((FP)1 + p.norm2());
    }

    GammaType getGamma() const { return static_cast<GammaType>(1.0) / invGamma; }

    MassType getMass() const { return mass; }
    void setMass(MassType newMass) { mass = newMass; }
    
    ChargeType getCharge() const { return charge; }
    void setCharge(ChargeType newCharge) { charge = newCharge; }

    FactorType getFactor() const { return factor; }
    void setFactor(FactorType newFactor) { factor = newFactor; }

private:

    PositionType position;
    MomentumType p;
    MassType mass;
    ChargeType charge;
    FactorType factor;
    GammaType invGamma;
};

typedef Particle<One> Particle1d;
typedef Particle<Two> Particle2d;
typedef Particle<Three> Particle3d;


} // namespace pica


#endif
