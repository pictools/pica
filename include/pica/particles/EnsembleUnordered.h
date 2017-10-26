#ifndef PICA_ENSEMBLEUNORDERED_H
#define PICA_ENSEMBLEUNORDERED_H


#include "pica/particles/ParticleTraits.h"


namespace pica {
    

// Representation of particles in an unordered array
// Particle are stored by ParticleArray class
template<class ParticleArray>
class EnsembleUnordered {
public:
    typedef typename ParticleArray::Particle Particle;
    typedef typename ParticleArray::ParticleRef ParticleRef;
    typedef typename ParticleArray::ConstParticleRef ConstParticleRef;
    typedef typename ParticleTraits<Particle>::PositionType PositionType;

    EnsembleUnordered(PositionType minPosition, PositionType maxPosition);

    PositionType getMinPosition() const;
    PositionType getMaxPosition() const;

    int size() const;

    template<class ConstParticleRef>
    void add(ConstParticleRef particle)
    {
        particles.pushBack(particle);
    }

    // Specific to this class
    ParticleArray& getParticles();
    const ParticleArray& getParticles() const;

protected:
    ParticleArray particles;
    PositionType minPosition, maxPosition;
};

template<class ParticleArray>
EnsembleUnordered<ParticleArray>::EnsembleUnordered(
        typename EnsembleUnordered<ParticleArray>::PositionType minPosition,
        typename EnsembleUnordered<ParticleArray>::PositionType maxPosition):
    minPosition(minPosition),
    maxPosition(maxPosition)
{
}

template<class ParticleArray>
typename EnsembleUnordered<ParticleArray>::PositionType EnsembleUnordered<ParticleArray>::getMinPosition() const
{ 
    return minPosition;
}

template<class ParticleArray>
typename EnsembleUnordered<ParticleArray>::PositionType EnsembleUnordered<ParticleArray>::getMaxPosition() const
{
    return maxPosition;
}

template<class ParticleArray>
int EnsembleUnordered<ParticleArray>::size() const
{ 
    return particles.size();
}

template<class ParticleArray>
ParticleArray& EnsembleUnordered<ParticleArray>::getParticles()
{ 
    return particles;
}

template<class ParticleArray>
const ParticleArray& EnsembleUnordered<ParticleArray>::getParticles() const
{ 
    return particles;
}


} // namespace pica


#endif
