#ifndef PICA_BORISPUSHER_H
#define PICA_BORISPUSHER_H


#include "pica/math/Constants.h"
#include "pica/math/Vectors.h"
#include "pica/particles/ParticleTraits.h"


namespace pica {


template<class Particle>
struct BorisPusher {

    typedef typename ParticleTraits<Particle>::PositionType PositionType;
    typedef typename ParticleTraits<Particle>::MomentumType MomentumType;

    template<class ParticleRef, typename Real>
    void push(ParticleRef* particle, const MomentumType& e, const MomentumType& b, Real dt)
    {
        const Real eCoeff = particle->getCharge() * dt /
            ((Real)2.0 * particle->getMass() * Constants<Real>::c());
        // The code below uses precomputed coefficient:
        // eCoeff = q * dt / (2 * m * c)
        MomentumType eMomentum = e * eCoeff;
        MomentumType um = particle->getP() + eMomentum;
        MomentumType t = b * eCoeff / sqrt((Real)1.0 + um.norm2());
        MomentumType uprime = um + cross(um, t);
        MomentumType s = t * (Real)2.0 / ((Real)1.0 + t.norm2());
        particle->setP((um + cross(uprime, s) + eMomentum));
        PositionType position = particle->getPosition();
        position += particle->getVelocity() * dt;
        particle->setPosition(position);
    }
};


} // namespace pica


#endif
