#ifndef PICA_PARTICLEARRAY_H
#define PICA_PARTICLEARRAY_H


#include "pica/math/Dimension.h"
#include "pica/particles/Particle.h"

#include <vector>


namespace pica {


    
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
template<typename PositionType, typename MomentumType>
class ParticleArraySoA {
public:

    class ParticleRef {
    public:
        ParticleRef(ParticleArraySoA& particleArray, int idx) :
            particleArray(particleArray),
            idx(idx)
        {}

        // ... (Implementation of particle-compatible interface)
    private:
        ParticleArraySoA& particleArray;
        int idx;
    };
    typedef const ParticleRef ConstParticleRef;

    int size() const { return positions.size(); }

    ParticleRef operator[](int idx) { return ParticleRef(this, idx); }
    ConstParticleRef operator[](int idx) const { return ConstParticleRef(this, idx); }

    void pushBack(ConstParticleRef particle) { /*particles.push_back(particle);*/ throw; }

private:

    std::vector<ScalarType<PositionType> > positions[VectorDimensionHelper<PositionType>::dimension];
    std::vector<ScalarType<MomentumType> > momentums[VectorDimensionHelper<MomentumType>::dimension];

    friend class ParticleRef;
    friend class ConstParticleRef;
};


} // namespace pica


#endif
