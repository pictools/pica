#ifndef PICA_PARTICLEARRAY_H
#define PICA_PARTICLEARRAY_H


#include "pica/math/Dimension.h"
#include "pica/math/Vectors.h"
#include "pica/particles/Particle.h"
#include "pica/particles/ParticleTraits.h"

#include <map>
#include <vector>
#include <string>


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
        typedef typename ParticleTraits<ParticleType>::FactorType FactorType;
        typedef typename ParticleTraits<ParticleType>::TypeIndexType TypeIndexType;
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
                result[d] = particles.ps[d][idx];
            return result * Constants<GammaType>::c() * pica::ParticleTypes::types[particles.typeIndex[idx]].mass;
        }

        MomentumType getP() const {
            MomentumType result;
            for (int d = 0; d < momentumDimension; d++)
                result[d] = particles.ps[d][idx];
            return result;
        }

        MomentumType getVelocity() const { return this->getP() * Constants<GammaType>::c() * particles.invGammas[idx]; }

        GammaType getGamma() const { return static_cast<GammaType>(1.0) / particles.invGammas[idx]; }

        MassType getMass() const { return ParticleTypes::types[particles.typeIndex[idx]].mass; }

        ChargeType getCharge() const { return ParticleTypes::types[particles.typeIndex[idx]].charge; }

        FactorType getFactor() const { return particles.factors[idx]; }

        TypeIndexType getType() const { return particles.typeIndex[idx]; }


    private:
        const ParticleArraySoA& particles;
        int idx;
    };

    class ParticleRef : public ConstParticleRef {
    public:

        typedef typename ParticleTraits<ParticleType>::PositionType PositionType;
        typedef typename ParticleTraits<ParticleType>::MomentumType MomentumType;
        typedef typename ParticleTraits<ParticleType>::GammaType GammaType;
        typedef typename ParticleTraits<ParticleType>::FactorType FactorType;
        typedef typename ParticleTraits<ParticleType>::TypeIndexType TypeIndexType;
        static const int dimension = VectorDimensionHelper<PositionType>::dimension;
        static const int momentumDimension = VectorDimensionHelper<MomentumType>::dimension;

        ParticleRef(ParticleArraySoA& particles, int idx) :
            ConstParticleRef(particles, idx),
            particles(particles),
            idx(idx)
        {}

        void setPosition(const PositionType& newPosition) {
            for (int d = 0; d < dimension; d++)
                particles.positions[d][idx] = newPosition[d];
        }

        void setMomentum(const MomentumType& newMomentum)
        {
            MomentumType p = newMomentum / (Constants<GammaType>::c() * ParticleTypes::types[particles.typeIndex[idx]].mass);
            for (int d = 0; d < momentumDimension; d++)
                particles.ps[d][idx] = p[d];
            particles.invGammas[idx] = static_cast<GammaType>(1.0) / sqrt(static_cast<GammaType>(1.0) + p.norm2());
        }

        void setP(const MomentumType& newP)
        {
            MomentumType p = newP;
            for (int d = 0; d < momentumDimension; d++)
                particles.ps[d][idx] = p[d];
            particles.invGammas[idx] = static_cast<GammaType>(1.0) / sqrt(static_cast<GammaType>(1.0) + p.norm2());
        }

        void setVelocity(const MomentumType& newVelocity)
        {
            MomentumType p = newVelocity / sqrt(constants::c * constants::c - newVelocity.norm2());
            for (int d = 0; d < momentumDimension; d++)
                particles.ps[d][idx] = p[d];
            particles.invGammas[idx] = static_cast<GammaType>(1.0) / sqrt(static_cast<GammaType>(1.0) + p.norm2());
        }

        void setFactor(FactorType newFactor) { particles.factors[idx] = newFactor; }

        void setType(TypeIndexType newTypeIndex) { particles.typeIndex[idx] = newTypeIndex; }

    private:
        // These intentionally overshadow members of ConstParticleRef
        ParticleArraySoA& particles;
        int idx;
    };

    int size() const { return static_cast<int>(typeIndex.size()); }

    ParticleRef operator[](int idx) { return ParticleRef(*this, idx); }
    ConstParticleRef operator[](int idx) const { return ConstParticleRef(*this, idx); }

    ParticleRef back() { return (*this)[size() - 1]; }
    ConstParticleRef back() const { return (*this)[size() - 1]; }

    template<class ConstParticleRefType>
    void pushBack(ConstParticleRefType particle)
    {
        const typename ParticleRef::PositionType position = particle.getPosition();
        for (int d = 0; d < ParticleRef::dimension; d++)
            positions[d].push_back(position[d]);
        const typename ParticleRef::MomentumType p = particle.getP();
        for (int d = 0; d < ParticleRef::momentumDimension; d++)
            ps[d].push_back(p[d]);
        factors.push_back(particle.getFactor());
        invGammas.push_back(static_cast<typename ParticleRef::GammaType>(1.0) / particle.getGamma());
        typeIndex.push_back(particle.getType());
    }
    void popBack()
    {
        for (int d = 0; d < ParticleRef::dimension; d++)
            positions[d].pop_back();
        for (int d = 0; d < ParticleRef::momentumDimension; d++)
            ps[d].pop_back();
        factors.pop_back();
        invGammas.pop_back();
        typeIndex.pop_back();
    }

private:

    std::vector<typename ScalarType<typename ParticleRef::PositionType>::Type> positions[ParticleRef::dimension];
    std::vector<typename ScalarType<typename ParticleRef::MomentumType>::Type> ps[ParticleRef::momentumDimension];
    std::vector<typename ParticleRef::FactorType> factors;
    std::vector<typename ParticleRef::GammaType> invGammas;
    std::vector<typename ParticleRef::TypeIndexType> typeIndex;

    friend class ParticleRef;
    friend class ConstParticleRef;
};


enum ParticleRepresentation { ParticleRepresentation_AoS, ParticleRepresentation_SoA };

inline std::string toString(ParticleRepresentation particleRepresentation)
{
    std::map<ParticleRepresentation, std::string> names;
    names[ParticleRepresentation_AoS] = "AoS";
    names[ParticleRepresentation_SoA] = "SoA";
    return names[particleRepresentation];
}


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
