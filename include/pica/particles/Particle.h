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
        charge(0.0) {}

    Particle(const PositionType& position, const MomentumType& momentum,
        MassType mass, ChargeType charge, FactorType factor = 1) :
        position(position), momentum(momentum), mass(mass), charge(charge), factor(factor) {}

    PositionType getPosition() const { return position; }
    void setPosition(const PositionType& newPosition) { position = newPosition; }

    MomentumType getMomentum() const { return momentum; }
    void setMomentum(const MomentumType& newMomentum) { momentum = newMomentum; }

    MomentumType getVelocity() const { return momentum / sqrt(sqr(getMass()) + (momentum / Constants<FP>::c()).norm2()); }
    void setVelocity(const MomentumType& newVelocity) { momentum = newVelocity * getMass() / sqrt(static_cast<Real>(1.0) - (newVelocity / constants::c).norm2()); }

    GammaType getGamma() const { return sqrt(static_cast<FP>(1.0) + (momentum / (getMass() * Constants<GammaType>::c())).norm2()); }

    MassType getMass() const { return mass; }
    void setMass(MassType newMass) { mass = newMass; }
    
    ChargeType getCharge() const { return charge; }
    void setCharge(ChargeType newCharge) { charge = newCharge; }

    FactorType getFactor() const { return factor; }
    void setFactor(FactorType newFactor) { factor = newFactor; }

private:

    PositionType position;
    MomentumType momentum;
    MassType mass;
    ChargeType charge;
    FactorType factor;

};

typedef Particle<One> Particle1d;
typedef Particle<Two> Particle2d;
typedef Particle<Three> Particle3d;


} // namespace pica


#endif
