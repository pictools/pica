#ifndef PICA_PARTICLE_INTERFACE_H
#define PICA_PARTICLE_INTERFACE_H

// This file serves as a documentation and should not be included as it will fail to build

namespace pica {


// This class defines an interface for particle classes to follow
class ParticleInterface {
public:

    // Either provide the following types or specialization of ParticleTraits
    typedef Vector3<double> PositionType; // or Vector2<double> or double
    typedef Vector3<double> MomentumType;
    typedef double GammaType;
    typedef double MassType;
    typedef double ChargeType;
    typedef double FactorType;

    ParticleInterface();
    ParticleInterface(const PositionType& position, const MomentumType& momentum,
        MassType mass, ChargeType charge, FactorType factor);

    MassType getMass() const;

    ChargeType getCharge() const;

    PositionType getPosition() const;
    void setPosition(const PositionType& newPosition);

    MomentumType getMomentum() const;
    void setMomentum(const MomentumType& newMomentum);

    MomentumType getVelocity() const;
    void setVelocity(const MomentumType& newVelocity);

    GammaType getGamma() const;
    
    FactorType getFactor() const;
    void setFactor(FactorType newFactor);
};


} // namespace pica


#endif

