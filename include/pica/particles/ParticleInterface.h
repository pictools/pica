#ifndef PICA_PARTICLE_INTERFACE_H
#define PICA_PARTICLE_INTERFACE_H


namespace pica {


// This class defines an interface for particle classes to follow
class ParticleInterface {
public:

    // Either provide the following types or specialization of ParticleTraits
    typedef FP3 PositionType;
    typedef FP3 MomemtumType;
    typedef FP GammaType;
    typedef FP ChargeType;
    typedef FP MassType;
    typedef FP FactorType;
    typedef int SpeciesType;

    ParticleInterface();
    ParticleInterface(const PositionType& position, const MomemtumType& momentum,
        const SpeciesType& species, FactorType factor);

    MassType getMass() const;
    ChargeType getCharge() const;
    
    SpeciesType getSpecies() const;
    void setSpecies(const SpeciesType& newSpecies);

    PositionType getPosition() const;
    void setPosition(const PositionType& newPosition);

    MomemtumType getMomentum() const;
    void setMomentum(const MomemtumType& newMomentum);

    MomemtumType getVelocity() const;
    void setVelocity(const MomemtumType& newVelocity);

    GammaType getGamma() const;
    
    FactorType getFactor() const;
    void setFactor(FactorType newFactor);
};


} // namespace pica


#endif

