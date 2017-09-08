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
class Particle_ {
public:

    // Types for conforming ParticleInterface
    typedef typename VectorTypeHelper<dimension, FP>::Type PositionType;
    typedef FP3 MomentumType;
    typedef FP GammaType;
    typedef FP ChargeType;
    typedef FP MassType;
    typedef FP FactorType;

    class SpeciesType {
    public:
        MassType getMass() const { return mass; }
        ChargeType getCharge() const { return charge; }
        std::string getName() const { return name; }
    private:
        ChargeType charge;
        MassType mass;
        std::string name;
    };

    Particle_() :
        factor(1) {}

    Particle_(const PositionType& position, const MomentumType& momentum,
        const SpeciesType& species, FactorType factor) :
        position(position), momentum(momentum), species(species), factor(factor) {}

    MassType getMass() const { return getSpecies().getMass(); }
    ChargeType getCharge() const { return getSpecies().getCharge(); }

    SpeciesType getSpecies() const { return species; }
    void setSpecies(const SpeciesType& newSpecies) { species = newSpecies; }

    PositionType getPosition() const { return position; }
    void setPosition(const PositionType& newPosition) { position = newPosition; }

    MomentumType getMomentum() const { return momentum; }
    void setMomentum(const MomentumType& newMomentum) { momentum = newMomentum; }

    MomentumType getVelocity() const { return momentum / sqrt(sqr(getMass()) + (momentum / Constants<FP>::c()).norm2()); }
    void setVelocity(const MomentumType& newVelocity) { momentum = getMass() * Constants<FP>::c() * newVelocity / sqrt(sqr(constants::c * constants::c) - newVelocity.norm2()); }

    GammaType getGamma() const { return sqrt(static_cast<FP>(1) / (static_cast<FP>(1) + (momentum / (getMass() * constants::c)).norm2())); }

    FactorType getFactor() const { return factor; }
    void setFactor(FactorType newFactor) { factor = newFactor; }

private:

    PositionType position;
    MomentumType momentum;
    SpeciesType species;
    FactorType factor;

};

typedef Particle_<One> Particle1d;
typedef Particle_<Two> Particle2d;
typedef Particle_<Three> Particle3d;

struct ParticleType
{
    FP mass;
    FP charge;

    // MDK-compatible
    std::string _name;
    std::string& name() { return _name; }
    const std::string& name() const { return _name; }

};


// Particle is internal particle representation.
class Particle
{
public:

    Particle()
    {
        typeIndex = 0;
        invGamma = 1;
        factor = 1;
        id = 0;
    }

    Particle(const FP3 & position, const FP3 & momentum, int _typeIndex, float _factor = 1.0f)
    {
        setType(_typeIndex);
        setFactor(_factor);
        setPosition(position);
        setMomentum(momentum);
        id = 0;
    }

    FP mass() const
    {
        return types[typeIndex].mass;
    }

    FP charge() const
    {
        return types[typeIndex].charge;
    }

    int getType() const
    {
        return typeIndex;
    }

    void setType(int newType)
    {
        typeIndex = newType;
    }

    const FP3 getPosition() const
    {
        return coords;
    }

    void setPosition(const FP3 & newPosition)
    {
        coords = newPosition;
    }

    const FP3 getMomentum() const
    {
        return p * Constants<FP>::c() * mass();
    }

    void setMomentum(const FP3 & newMomentum)
    {
        p = newMomentum / (Constants<FP>::c() * mass());
        invGamma = (FP)1 / sqrt((FP)1 + p.norm2());
    }

    const FP3 getVelocity() const
    {
        return p * (Constants<FP>::c() * invGamma);
    }

    void setVelocity(const FP3 & newVelocity)
    {
        p = newVelocity / sqrt(Constants<FP>::c() * Constants<FP>::c() - newVelocity.norm2());
        invGamma = (FP)1 / sqrt((FP)1 + p.norm2());
    }

    FP gamma() const
    {
        return (FP)1 / invGamma;
    }

    FP getFactor() const
    {
        return factor;
    }

    void setFactor(FP newFactor)
    {
        factor = (float)newFactor;
    }

    static void getOffsets(size_t* offsets);

    static std::vector<ParticleType> typesVector;
    static const ParticleType * types; // raw pointer to elements of typesVector
    static int numTypes; // size of typesVector

public:

    FP3 coords;
    FP invGamma; // invGamma = 1/sqrt(1+p^2), v = p * c * invGamma
    FP3 p; // physical momentum divided by (m * c)
    short typeIndex;
    short id;
    float factor;

    friend class CurrentDeposition;
    friend class ParticlePush;
    friend class ParticleTracking;
    friend class Checkpoint;
};

} // namespace pica


#endif
