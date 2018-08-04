#ifndef PICA_BORISPUSHER_H
#define PICA_BORISPUSHER_H


#include "pica/math/Constants.h"
#include "pica/math/Vectors.h"
#include "pica/particles/ParticleTraits.h"


namespace pica {


template<class Particle, typename Real>
struct BorisPusher {

    typedef typename ParticleTraits<Particle>::PositionType PositionType;
    typedef typename ParticleTraits<Particle>::MomentumType MomentumType;

    std::vector<Real> coeffVector;
    Real* coeff;
    Real dt;

    BorisPusher(Real dt)
        : dt(dt)
    {
        coeffVector.resize(ParticleTypes::numTypes);
        coeff = &coeffVector[0];
        for (int i = 0; i < ParticleTypes::numTypes; i++)
            coeff[i] = ParticleTypes::typesVector[i].charge * dt / (2.0 * ParticleTypes::typesVector[i].mass * Constants<FP>::c());
    }

    template<class ParticleRef>
    void push(ParticleRef* particle, const MomentumType& e, const MomentumType& b)
    {
        // The code below uses precomputed coefficient:
        // eCoeff = q * dt / (2 * m * c)
        Real eCoeff = coeff[particle->getType()];
        MomentumType eMomentum = e * eCoeff;
        MomentumType um = particle->getP() + eMomentum;
        MomentumType t = b * (eCoeff / sqrt((Real)1.0 + um.norm2()));
        MomentumType uprime = um + cross(um, t);
        MomentumType s = t * ((Real)2.0 / ((Real)1.0 + t.norm2()));
        particle->setP(um + cross(uprime, s) + eMomentum);
        PositionType position = particle->getPosition();
        position += particle->getVelocity() * dt;
        particle->setPosition(position);
    }
};


} // namespace pica


#endif
