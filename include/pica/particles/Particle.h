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

typedef FP Real;
typedef Real MassType;
typedef Real ChargeType;

struct ParticleType {
    MassType mass;
    ChargeType charge;
};

namespace ParticleTypes
{
    extern std::vector<ParticleType> typesVector;
    extern const ParticleType* types;
    extern int numTypes;
};

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
    typedef Real FactorType;
    typedef short TypeIndexType;

    Particle() :
        factor(1),
        invGamma(1.0),
        typeIndex(0)
    {}

    Particle(const PositionType& position, const MomentumType& momentum,
        FactorType factor = 1, TypeIndexType typeIndex = 0) :
        position(position), factor(factor), typeIndex(typeIndex)
    {
        setMomentum(momentum);
    }

    PositionType getPosition() const { return position; }
    void setPosition(const PositionType& newPosition) { position = newPosition; }

    MomentumType getMomentum() const
    {
        return p * Constants<MassType>::c() * getMass();
    }

    void setMomentum(const MomentumType& newMomentum)
    {
        p = newMomentum / (Constants<GammaType>::c() * getMass());
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

    MassType getMass() const { return ParticleTypes::types[typeIndex].mass; }

    ChargeType getCharge() const { return ParticleTypes::types[typeIndex].charge; }

    FactorType getFactor() const { return factor; }
    void setFactor(FactorType newFactor) { factor = newFactor; }

    TypeIndexType getType() const { return typeIndex; }
    void setType(TypeIndexType newType) { typeIndex = newType; }

private:

    PositionType position;
    MomentumType p;
    FactorType factor;
    GammaType invGamma;
    TypeIndexType typeIndex;
};

typedef Particle<One> Particle1d;
typedef Particle<Two> Particle2d;
typedef Particle<Three> Particle3d;


} // namespace pica


#endif
