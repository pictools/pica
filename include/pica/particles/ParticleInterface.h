#ifndef PICA_PARTICLE_INTERFACE_H
#define PICA_PARTICLE_INTERFACE_H

// This file serves as a documentation and should not be included as it will fail to build

namespace pica {


// This class defines an interface for particle classes to follow
class ParticleInterface {
public:

    // Either provide the following types or specialization of ParticleTraits
    typedef Vector3<double> PositionType;
    typedef Vector3<double> MomentumType;
    typedef double GammaType;
    typedef double ChargeType;
    typedef double MassType;
    typedef double FactorType;
    typedef int SpeciesType;

    ParticleInterface();
    ParticleInterface(const PositionType& position, const MomentumType& momentum,
        const SpeciesType& species, FactorType factor);

    MassType getMass() const;
    ChargeType getCharge() const;
    
    SpeciesType getSpecies() const;
    void setSpecies(const SpeciesType& newSpecies);

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

