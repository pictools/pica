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

struct ParticleType
{
    FP mass;
    FP charge;

    // MDK-compatible
    std::string _name;
    std::string& name() { return _name; }
    const std::string& name() const { return _name; }

    // Ionization
    int numChargeStates;
    FP chargeStep;
    int ionizationL;
    FP ionizationPotential;
    int numElectrons;
    int z0;

    ParticleType(): numChargeStates(1), chargeStep(1), ionizationL(0),
        ionizationPotential(0), numElectrons(0), z0(0) {}
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
        return p * constants::c * mass();
    }

    void setMomentum(const FP3 & newMomentum)
    {
        p = newMomentum / (constants::c * mass());
        invGamma = (FP)1 / sqrt((FP)1 + p.norm2());
    }

    const FP3 getVelocity() const
    {
        return p * (constants::c * invGamma);
    }

    void setVelocity(const FP3 & newVelocity)
    {
        p = newVelocity / sqrt(constants::c * constants::c - newVelocity.norm2());
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
