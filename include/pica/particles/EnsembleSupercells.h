#ifndef PICA_ENSEMBLESUPERCELLS_H
#define PICA_ENSEMBLESUPERCELLS_H


#include "pica/math/Vectors.h"
#include "pica/particles/ParticleTraits.h"
#include "pica/utility/Array.h"


namespace pica {


// Representation of particles in an unordered array
// Particle are stored by ParticleArray class
template<class ParticleArray>
class EnsembleSupercells {
public:
    typedef typename ParticleArray::Particle Particle;
    typedef typename ParticleArray::ParticleRef ParticleRef;
    typedef typename ParticleArray::ConstParticleRef ConstParticleRef;
    typedef typename ParticleTraits<Particle>::PositionType PositionType;
    static const int dimension = VectorDimensionHelper<PositionType>::dimension;
    typedef typename VectorIntTypeHelper<dimension, int>::Type IndexType;

    EnsembleSupercells(PositionType minPosition, PositionType maxPosition,
        IndexType numCells, IndexType numCellsPerSupercell) :
        minPosition(minPosition), maxPosition(maxPosition), numCells(numCells),
        numCellsPerSupercell(numCellsPerSupercell) {}

    PositionType getMinPosition() const { return minPosition; }
    PositionType getMaxPosition() const { return maxPosition; }

    int size() const { return 0; /* !!! */ }

    template<class ConstParticleRef>
    void add(ConstParticleRef particle)
    {
        /* !!! */
    }

    // Specific to this class
    IndexType getNumSupercells() const { return supercells.getSize(); }
    ParticleArray& getParticles(IndexType supercellIdx) { return supercells(supercellIdx); }
    const ParticleArray& getParticles(IndexType supercellIdx) const { return supercells(supercellIdx); }

private:
    PositionType minPosition;
    PositionType maxPosition;
    IndexType numCells;
    IndexType numCellsPerSupercell;

    typedef typename ArrayIntTypeHelper<dimension, ParticleArray>::Type Supercells;
    Supercells supercells;

};


} // namespace pica


#endif
