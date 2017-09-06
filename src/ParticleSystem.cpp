#include "pica/ParticleSystem.h"

#include "pica/utility/Utility.h"


namespace pica {

ParticleSystem::ParticleSystem(int allocatedParticles)
{
    particles.reserve(allocatedParticles);
    begin();
}


ParticleSystem::ParticleSystem(const ParticleSystem & system)
{
    particles = system.particles;
    currentParticle = raw() + (system.currentParticle - system.raw());
}


ParticleSystem& ParticleSystem::operator =(const ParticleSystem & system)
{
    particles = system.particles;
    currentParticle = raw() + (system.currentParticle - system.raw());
    return *this;
}


Particle* ParticleSystem::raw()
{
    return ptr(particles);
}


const Particle* ParticleSystem::raw() const
{
    return ptr(particles);
}


void ParticleSystem::add(const Particle & newParticle)
{
    Particle * oldRaw = raw();
    particles.push_back(newParticle);
    currentParticle = raw() + (currentParticle - oldRaw);
}

} // namespace pica
