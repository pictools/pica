#ifndef PICA_PARTICLEARRAY_H
#define PICA_PARTICLEARRAY_H


#include "pica/math/Dimension.h"
#include "pica/particles/Particle.h"
#include "pica/particles/ParticleTraits.h"

#include <vector>


namespace pica {


enum ParticleRepresentation { ParticleRepresentation_AoS, ParticleRepresentation_SoA };

// Traits class to provide a Type corresponding to array of particles
// according to the given representation
template<class Particle, ParticleRepresentation storage>
struct ParticleArray {
};

template<class Particle, ParticleRepresentation_AoS>
struct ParticleArray {
    typedef ParticleArrayAoS<Particle> Type;
};

template<class Particle, ParticleRepresentation_SoA>
struct ParticleArray {
    typedef ParticleArraySoA<Particle> Type;
};

  
// Collection of particles with array-like semantics,
// representation as array of structures
template<class Particle>
class ParticleArrayAoS {
public:
    typedef Particle& ParticleRef;
    typedef const ParticleRef ConstParticleRef;

    int size() const { return particles.size(); }

    ParticleRef operator[](int idx) { return particles[idx]; }
    ConstParticleRef operator[](int idx) const { return particles[idx]; }

    void pushBack(ConstParticleRef particle) { particles.push_back(particle); }

private:

    std::vector<Particle> particles;
};


// Collection of particles with array-like semantics,
// representation as structure of arrays
template<class Particle>
class ParticleArraySoA {
public:

    class ConstParticleRef {
    public:
        ConstParticleRef(const ParticleArraySoA& particles, int idx):
            particles(particles),
            idx(idx)
        {}

        typedef typename ParticleTraits<Particle>::PositionType PositionType;
        typedef typename ParticleTraits<Particle>::PositionType MomentumType;
        typedef typename ParticleTraits<Particle>::GammaType GammaType;
        typedef typename ParticleTraits<Particle>::MassType MassType;
        typedef typename ParticleTraits<Particle>::ChargeType ChargeType;
        typedef typename ParticleTraits<Particle>::FactorType FactorType;
        const int dimension = VectorDimensionHelper<PositionType>::dimension;
        const int momentumDimension = VectorDimensionHelper<MomentumType>::dimension;

        PositionType getPosition() const {
            PositionType result;
            for (int d = 0; d < dimension; d++)
                result[d] = particles.positions[d][idx];
            return result;
        }

        MomentumType getMomentum() const {
            MomentumType result;
            for (int d = 0; d < momentumDimension; d++)
                result[d] = particles.momentums[d][idx];
            return result;
        }

        MomentumType getVelocity() const { return getMomentum() / sqrt(sqr(getMass()) + (momentum / Constants<FP>::c()).norm2()); }
 
        GammaType getGamma() const { return sqrt(static_cast<GammaType>(1) / (static_cast<GammaType>(1) + (getMomentum() / (getMass() * constants::c)).norm2())); }

        MassType getMass() const { return particles.masses[idx]; }
 
        ChargeType getCharge() const { return particles.charges[idx]; }
 
        FactorType getFactor() const { return particles.factors[idx]; }

    private:
        const ParticleArraySoA& particles;
        int idx;
    };

    class ParticleRef : public ConstParticleRef {
    public:

        ParticleRef(ParticleArraySoA& particles, int idx) :
            ConstParticleRef(particles, idx),
            particles(particles),
            idx(idx)
        {}

        void setPosition(const PositionType& newPosition) {
            for (int d = 0; d < dimension; d++)
                particles.positions[d][idx] = newPosition[d];
        }

        void setMomentum(const MomentumType& newMomentum) { 
            for (int d = 0; d < momentumDimension; d++)
                particles.momentums[d][idx] = newMomentum[d];
        }

        void setVelocity(const MomentumType& newVelocity) { setMomentum(getMass() * Constants<GammaType>::c() * newVelocity / sqrt(sqr(constants::c * constants::c) - newVelocity.norm2())); }

        void setMass(MassType newMass) { particles.masses[idx] = newMass; }

        void setCharge(ChargeType newCharge) { particles.charges[idx] = newCharge; }

        void setFactor(FactorType newFactor) { particles.factors[idx] = newFactor; }

    private:
        // These intentionally overshadow members of ConstParticleRef
        ParticleArraySoA& particles;
        int idx;
    };

    int size() const { return positions.size(); }

    ParticleRef operator[](int idx) { return ParticleRef(this, idx); }
    ConstParticleRef operator[](int idx) const { return ConstParticleRef(this, idx); }

    template<class ConstParticleRefType>
    void pushBack(ConstParticleRefType particle)
    {
        const ParticleRef::PositionType position = particle.getPosition();
        for (int d = 0; d < dimension; d++)
            positions[d].push_back = position[d];
        const ParticleRef::MomentumType momentum = particle.getMomentum();
        for (int d = 0; d < momentumDimension; d++)
            momentums[d].push_back = momentum[d];
        masses.push_back(particle.getMass());
        charges.push_back(particle.getCharge());
        factors.push_back(particle.getFactor());
    }

private:

    std::vector<ScalarType<ParticleRef::PositionType> > positions[ParticleRef::dimension];
    std::vector<ScalarType<ParticleRef::MomentumType> > momentums[ParticleRef::momentumDimension];
    std::vector<ParticleRef::MassType> masses;
    std::vector<ParticleRef::ChargeType> charges;
    std::vector<ParticleRef::FactorType> factors;

    friend class ParticleRef;
    friend class ConstParticleRef;
};


} // namespace pica


#endif
