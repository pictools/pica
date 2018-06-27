#ifndef PICA_PARTICLE_BASELINE_H
#define PICA_PARTICLE_BASELINE_H


#include "pica/math/Constants.h"
#include "pica/math/FP.h"
#include "pica/math/Vectors.h"

#include <cmath>
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

// Baseline representation for particle
template<Dimension dimension>
class ParticleBaseline {
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

    ParticleBaseline() :
        factor(1),
        typeIndex(0)
    {}

    ParticleBaseline(const PositionType& position, const MomentumType& momentum, FactorType factor = 1, TypeIndexType typeIndex = 0) :
        position(position), momentum(momentum), mass(mass), charge(charge), factor(factor), typeIndex(typeIndex) {}

    PositionType getPosition() const { return position; }
    void setPosition(const PositionType& newPosition) { position = newPosition; }

    MomentumType getMomentum() const { return momentum; }
    void setMomentum(const MomentumType& newMomentum) { momentum = newMomentum; }

    MomentumType getVelocity() const { return momentum / sqrt(sqr(getMass()) + (momentum / Constants<FP>::c()).norm2()); }
    void setVelocity(const MomentumType& newVelocity) { momentum = newVelocity * getMass() / sqrt(static_cast<Real>(1.0) - (newVelocity / constants::c).norm2()); }

    GammaType getGamma() const { return sqrt(static_cast<FP>(1.0) + (momentum / (getMass() * Constants<GammaType>::c())).norm2()); }

    MassType getMass() const { return ParticleTypes::types[typeIndex].mass; }

    ChargeType getCharge() const { return ParticleTypes::types[typeIndex].charge; }

    FactorType getFactor() const { return factor; }
    void setFactor(FactorType newFactor) { factor = newFactor; }

    TypeIndexType getType() const { return typeIndex; }
    void setType(TypeIndexType newType) { typeIndex = newType; }

private:

    PositionType position;
    MomentumType momentum;
    FactorType factor;
    TypeIndexType typeIndex;
};


} // namespace pica


#endif