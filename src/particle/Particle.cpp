#include "pica/particle/Particle.h"


namespace pica {

std::vector<ParticleType> Particle::typesVector;
const ParticleType * Particle::types;
int Particle::numTypes;


void Particle::getOffsets(size_t* offsets)
{
    offsets[0] = offsetof(Particle, coords);
    offsets[1] = offsetof(Particle, invGamma);
    offsets[2] = offsetof(Particle, p);
    offsets[3] = offsetof(Particle, typeIndex);
    offsets[4] = offsetof(Particle, id);
    offsets[5] = offsetof(Particle, factor);
}

} // namespace pica
