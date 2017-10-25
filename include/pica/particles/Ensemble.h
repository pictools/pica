#ifndef PICA_ENSEMBLE_H
#define PICA_ENSEMBLE_H


#include "pica/particles/Particle.h"
#include "pica/particles/ParticleTraits.h"
#include "pica/particles/ParticleSystem.h"
#include "pica/threading/OpenMPHelper.h"

#include <algorithm>
#include <vector>

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

    EnsembleUnordered(PositionType minPosition, PositionType maxPosition) :
        minPosition(minPosition),
        maxPosition(maxPosition)
    {}

    int size() const { return particles.size(); }

    template<class ConstParticleRef>
    void add(ConstParticleRef particle) { particles.pushBack(particle); }

    // Specific to this class
    ParticleRef operator[](int idx) { return particles[idx]; }
    ConstParticleRef operator[](int idx) const { return particles[idx]; }

protected:
    ParticleArray particles;
    PositionType minPosition, maxPosition;
};


// Representation of particles in an unordered array
// Particle are stored by ParticleArray class
template<class ParticleArray>
class EnsembleOrdered : public EnsembleUnordered<ParticleArray> {
public:
    EnsembleOrdered(PositionType minPosition, PositionType maxPosition):
        EnsembleUnordered(minPosition, maxPosition)
    {}

    void reorder()
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
};


enum EnsembleRepresentation { EnsembleRepresentation_Unordered, EnsembleRepresentation_Ordered };

// Traits class to provide a Type corresponding to array of particles
// according to the given representation
template<class ParticleArray, EnsembleRepresentation storage>
struct Ensemble_ {
};

template<class ParticleArray>
struct Ensemble_<ParticleArray, EnsembleRepresentation_Unordered> {
    typedef EnsembleUnordered<ParticleArray> Type;
};

template<class ParticleArray>
struct Ensemble_<ParticleArray, EnsembleRepresentation_Ordered> {
    typedef EnsembleOrdered<ParticleArray> Type;
};


// Class for accessing to particles via internal iterator.
class Ensemble
{
public:

    Particle * currentParticle; // Pointer to particle where iterator currently is
    std::vector<ParticleType> &type;  // MDK-compatible vector of types

    Ensemble(const Int3& numInternalCells, const FP3& minCoords, const FP3& steps, int dimensionality);

    Ensemble(const Ensemble &ensemble);
    Ensemble& operator =(const Ensemble &ensemble);

    ~Ensemble();

    void add(const Particle& newParticle);
    void addParticle(const Particle& newParticle) { add(newParticle); }

    bool isMigratingParticle(const Particle& p, const Int3& idx) const
    { return (getSystemIdx(p) != idx); }

    void begin();
    void next();
    bool end();

    void deleteCurrentParticle();// Delete current particle and make currentParticle invalid
                                 // until call to next() that will set it to the next particle

    ParticleSystem * system(int i, int j, int k)
    {
        return &systems[k + numSystems.z * (j + numSystems.y * i)];
    }

    ParticleSystem * system(const Int3 & idx)
    {
        return system(idx.x, idx.y, idx.z);
    }

    const ParticleSystem * system(int i, int j, int k) const
    {
        return &systems[k + numSystems.z * (j + numSystems.y * i)];
    }

    const ParticleSystem * system(const Int3 & idx) const
    {
        return system(idx.x, idx.y, idx.z);
    }

    Int3 numSystems;

    void clear();

    int getNumAllocatedParticles() const;
    int getNumParticles() const;
    void shrinkToFit(); // decrease amount of memory allocated for particles, non-imperative (same as in std::vector) 
    void tryShrink();

    // MDK-compatible
    Int3 getCellIndex(const Particle& p) const
    { return getSystemIdx(p); }

    void dumpParticles(Particle* particles, const Int3 * minCellIdx,
        const Int3 * maxCellIdx) const;
    void dumpParticleInCellCount(int * particleInCellCount,
        const Int3 * minCellIdx, const Int3 * maxCellIdx) const;
    void dumpOutgoingParticles(std::vector<Particle>& particles);

    bool inParticlesLoop;
    void addPendingParticles();

     /* Get particle type index. First look for case-sensitive match,
     then for case-insensitive */
    int getTypeIdx(const std::string& name) const;

private:

    Int3 getSystemIdx(const Particle& p) const
    {  return truncate((p.getPosition() - origin) * invSteps); }

    ParticleSystem * pendingSystem(int i, int j, int k)
    {  return &pendingSystems[k + numSystems.z * (j + numSystems.y * i)]; }

    ParticleSystem * pendingSystem(Int3 idx)
    {  return pendingSystem(idx.x, idx.y, idx.z); }

    void computeBoundaryIndexes();

    std::vector<ParticleSystem> systems;
    std::vector<ParticleSystem> pendingSystems;
    std::vector<omp_lock_t> locks;
    omp_lock_t logMutex;
    int totalNumSystems;
    FP3 origin, invSteps;
    ParticleSystem * currentSystem;
    std::vector<std::vector<Particle> > pendingParticles;
    std::vector<int> boundaryCellIndexes;

    // temporary solution for 1D and 2D
    int dimensionality;
    FP3 centerPoint; // for 1D and 2D simulations
};

} // namespace pica


#endif
