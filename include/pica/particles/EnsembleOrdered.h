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
    EnsembleOrdered(PositionType minPosition, PositionType maxPosition);
    void reorder();
};

template<class ParticleArray>
EnsembleOrdered<ParticleArray>::EnsembleOrdered(
        typename EnsembleOrdered<ParticleArray>::PositionType minPosition,
        typename EnsembleOrdered<ParticleArray>::PositionType maxPosition):
    EnsembleUnordered<ParticleArray>(minPosition, maxPosition)
{
}

template<class ParticleArray>
void EnsembleOrdered<ParticleArray>::reorder()
{
    struct ParticleIndexComparator {
        ParticleIndexComparator(const ParticleArray& particles) :
            particles(particles),
            dimension(VectorDimensionHelper<PositionType>::dimension)
        {}

        bool operator()(int first, int second)
        {
            for (int d = 0; d < dimension; d++)
                if (particles[first].getPosition()[d] < particles[second].getPosition()[d])
                    return true;
                else if (particles[first].getPosition()[d] > particles[second].getPosition()[d])
                    return false;
            return false;
        }

    private:
        const ParticleArray& particles;
        const int dimension;
    };

    std::vector<int> indexes(particles.size());
    for (int i = 0; i < indexes.size(); i++)
        indexes[i] = i;
    std::sort(indexes.begin(), indexes.end(), ParticleIndexComparator(particles));
    ParticleArray newParticles;
    for (int i = 0; i < indexes.size(); i++)
        newParticles.pushBack(particles[indexes[i]]);
    particles = newParticles;
}


} // namespace pica


#endif
