#ifndef PICA_ENSEMBLE_H
#define PICA_ENSEMBLE_H


#include "pica/particles/Particle.h"
#include "pica/particles/ParticleSystem.h"
#include "pica/threading/OpenMPHelper.h"

#include <vector>

namespace pica {

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
