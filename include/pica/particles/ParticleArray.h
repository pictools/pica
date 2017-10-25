#ifndef PICA_PARTICLEARRAY_H
#define PICA_PARTICLEARRAY_H


#include "pica/math/Dimension.h"
#include "pica/math/Vectors.h"
#include "pica/particles/Particle.h"
#include "pica/particles/ParticleTraits.h"

#include <vector>


namespace pica {
   
  
// Collection of particles with array-like semantics,
// representation as array of structures
template<class ParticleType>
class ParticleArrayAoS {
public:
    typedef ParticleType Particle;
    typedef Particle& ParticleRef;
    typedef const Particle& ConstParticleRef;

    int size() const { return static_cast<int>(particles.size()); }

    ParticleRef operator[](int idx) { return particles[idx]; }
    ConstParticleRef operator[](int idx) const { return particles[idx]; }

    ParticleRef back() { return particles.back(); }
    ConstParticleRef back() const { return particles.back(); }

    void pushBack(ConstParticleRef particle) { particles.push_back(particle); }
    void popBack() { particles.pop_back(); }

private:

    std::vector<Particle> particles;
};


// Collection of particles with array-like semantics,
// representation as structure of arrays
template<class ParticleType>
class ParticleArraySoA {
public:

    typedef ParticleType Particle;

    class ConstParticleRef {
    public:
        ConstParticleRef(const ParticleArraySoA& particles, int idx):
            particles(particles),
            idx(idx)
        {}

        typedef typename ParticleTraits<ParticleType>::PositionType PositionType;
        typedef typename ParticleTraits<ParticleType>::MomentumType MomentumType;
        typedef typename ParticleTraits<ParticleType>::GammaType GammaType;
        typedef typename ParticleTraits<ParticleType>::MassType MassType;
        typedef typename ParticleTraits<ParticleType>::ChargeType ChargeType;
        typedef typename ParticleTraits<ParticleType>::FactorType FactorType;
        static const int dimension = VectorDimensionHelper<PositionType>::dimension;
        static const int momentumDimension = VectorDimensionHelper<MomentumType>::dimension;

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

        MomentumType getVelocity() const { return getMomentum() / sqrt(sqr(getMass()) + (getMomentum() / Constants<FP>::c()).norm2()); }
 
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

    int size() const { return static_cast<int>(masses.size()); }

    ParticleRef operator[](int idx) { return ParticleRef(*this, idx); }
    ConstParticleRef operator[](int idx) const { return ConstParticleRef(*this, idx); }

    ParticleRef back() { return (*this)[size() - 1]; }
    ConstParticleRef back() const { return (*this)[size() - 1]; }

    template<class ConstParticleRefType>
    void pushBack(ConstParticleRefType particle)
    {
        const ParticleRef::PositionType position = particle.getPosition();
        for (int d = 0; d < ParticleRef::dimension; d++)
            positions[d].push_back(position[d]);
        const ParticleRef::MomentumType momentum = particle.getMomentum();
        for (int d = 0; d < ParticleRef::momentumDimension; d++)
            momentums[d].push_back(momentum[d]);
        masses.push_back(particle.getMass());
        charges.push_back(particle.getCharge());
        factors.push_back(particle.getFactor());
    }
    void popBack()
    {
        for (int d = 0; d < ParticleRef::dimension; d++)
            positions[d].pop_back();
        for (int d = 0; d < ParticleRef::momentumDimension; d++)
            momentums[d].pop_back();
        masses.pop_back();
        charges.pop_back();
        factors.pop_back();
    }

private:

    std::vector<typename ScalarType<typename ParticleRef::PositionType>::Type> positions[ParticleRef::dimension];
    std::vector<typename ScalarType<typename ParticleRef::MomentumType>::Type> momentums[ParticleRef::momentumDimension];
    std::vector<typename ParticleRef::MassType> masses;
    std::vector<typename ParticleRef::ChargeType> charges;
    std::vector<typename ParticleRef::FactorType> factors;

    friend class ParticleRef;
    friend class ConstParticleRef;
};


enum ParticleRepresentation { ParticleRepresentation_AoS, ParticleRepresentation_SoA };

// Traits class to provide a Type corresponding to array of particles
// according to the given representation
template<class Particle, ParticleRepresentation storage>
struct ParticleArray {
};

template<class Particle>
struct ParticleArray<Particle, ParticleRepresentation_AoS> {
    typedef ParticleArrayAoS<Particle> Type;
};

template<class Particle>
struct ParticleArray<Particle, ParticleRepresentation_SoA> {
    typedef ParticleArraySoA<Particle> Type;
};



} // namespace pica


#endif
