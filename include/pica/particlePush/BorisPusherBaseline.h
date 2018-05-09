#ifndef PICA_BORISPUSHERBASELINE_H
#define PICA_BORISPUSHERBASELINE_H


#include "pica/math/Constants.h"
#include "pica/math/Vectors.h"
#include "pica/particles/ParticleTraits.h"


namespace pica {


template<class Particle>
struct BorisPusherBaseline {

    typedef typename ParticleTraits<Particle>::PositionType PositionType;
    typedef typename ParticleTraits<Particle>::MomentumType MomentumType;

    template<class ParticleRef, typename Real>
    void push(ParticleRef particle, const MomentumType& e, const MomentumType& b, Real dt)
    {
        MomentumType eMomentum = e * particle.getCharge() * dt /
            ((Real)2.0 * particle.getMass() * Constants<Real>::c());
        MomentumType um = particle.getMomentum() / (particle.getMass() * Constants<Real>::c()) + eMomentum;
        MomentumType t = b * particle.getCharge() * dt /
            ((Real)2.0 * particle.getMass() * Constants<Real>::c() * sqrt((Real)1.0 + um.norm2()));
        MomentumType uprime = um + cross(um, t);
        MomentumType s = t * (Real)2.0 / ((Real)1.0 + t.norm2());
        particle.setMomentum((um + cross(uprime, s) + eMomentum) * particle.getMass() * Constants<Real>::c());
        PositionType position = particle.getPosition();
        position += particle.getVelocity() * dt;
        particle.setPosition(position);
    }
};


} // namespace pica


#endif
