#ifndef PICA_PARTICLESYSTEM_H
#define PICA_PARTICLESYSTEM_H


#include "pica/Particle.h"

#include <vector>


namespace pica {

// System of particles with access as internal iterator provided
// via currentParticle, begin(), next(), end().
class ParticleSystem
{
public:

    /* Create a particle system with no particles and
    allocate memory for given number of particles */
    ParticleSystem(int allocatedSize = 1);
    ParticleSystem(const ParticleSystem & system);

    ParticleSystem& operator =(const ParticleSystem & system);

    void add(const Particle & newParticle);

    Particle * currentParticle;
    
    Particle * raw();
    const Particle * raw() const;

    void begin()
    {
        currentParticle = raw();
    }

    void next()
    {
        ++currentParticle;
    }

    bool end()
    {
        return ((currentParticle - raw()) >= size());
    }

    void deleteParticle(int i)
    {
        particles[i] = particles.back();
        particles.pop_back();
    }

    void deleteCurrentParticle()
    {
        *(currentParticle--) = particles.back();
        particles.pop_back();
    }

    void clear()
    {
        particles.clear();
    }

    int size() const
    {
        return (int)particles.size();
    }

private:

    std::vector<Particle> particles;

    friend class Ensemble;
    friend class CurrentDeposition;
    friend class Ionization;
    friend class ParticlePush;
};

} // namespace pica


#endif
