#ifndef PICA_ENSEMBLEORDERED_H
#define PICA_ENSEMBLEORDERED_H


#include "pica/particles/EnsembleUnordered.h"

#include <algorithm>
#include <vector>


namespace pica {


// Representation of particles in an unordered array
// Particle are stored by ParticleArray class
template<class ParticleArray>
class EnsembleOrdered : public EnsembleUnordered<ParticleArray> {
public:
    using typename EnsembleUnordered<ParticleArray>::PositionType;
    EnsembleOrdered(PositionType minPosition, PositionType maxPosition);
    void reorder();
    bool operator()(int first, int second) const;
private:
    int dimension;
};

template<class ParticleArray>
EnsembleOrdered<ParticleArray>::EnsembleOrdered(
        typename EnsembleOrdered<ParticleArray>::PositionType minPosition,
        typename EnsembleOrdered<ParticleArray>::PositionType maxPosition):
    EnsembleUnordered<ParticleArray>(minPosition, maxPosition)
{
    dimension = VectorDimensionHelper<PositionType>::dimension;
}

template<class ParticleArray>
bool EnsembleOrdered<ParticleArray>::operator()(int first, int second) const
{
    for (int d = 0; d < 3; d++)
        if (this->particles[first].getPosition()[d] < this->particles[second].getPosition()[d])
            return true;
        else if (this->particles[first].getPosition()[d] > this->particles[second].getPosition()[d])
            return false;
}


template<class ParticleArray>
void EnsembleOrdered<ParticleArray>::reorder()
{
    std::vector<int> indexes(this->particles.size());
    for (int i = 0; i < indexes.size(); i++)
        indexes[i] = i;
    std::sort(indexes.begin(), indexes.end(), *this);
    ParticleArray newParticles;
    for (int i = 0; i < indexes.size(); i++)
        newParticles.pushBack(this->particles[indexes[i]]);
    this->particles = newParticles;
}


} // namespace pica


#endif
