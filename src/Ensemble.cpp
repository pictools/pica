#include "pica/Ensemble.h"

#include "pica/OpenMPHelper.h"
#include "pica/Utility.h"

#include <algorithm>


namespace pica {

Ensemble::Ensemble(const Int3& numInternalCells, const FP3& minCoords, const FP3& steps, int _dimensionality):
    type(Particle::typesVector),
    numSystems(numInternalCells + Int3(2, 2, 2)),
    origin(minCoords - steps),
    dimensionality(_dimensionality)
{
    totalNumSystems = numSystems.x * numSystems.y * numSystems.z;
    // Add fake system in the end to correctly handle next() and end().
    systems.resize(totalNumSystems + 1);
    pendingSystems.resize(totalNumSystems + 1);
    locks.resize(totalNumSystems);
    for (int i = 0; i < totalNumSystems; i++)
        omp_init_lock(&locks[i]);
    omp_init_lock(&logMutex);
    invSteps = FP3(1.0, 1.0, 1.0) / steps;
    pendingParticles.resize(omp_get_max_threads());
    inParticlesLoop = false;
    computeBoundaryIndexes();
    centerPoint = minCoords + 0.5 * steps * numInternalCells;
}


Ensemble::Ensemble(const Ensemble &ensemble):
    systems(ensemble.systems),
    pendingSystems(ensemble.pendingSystems),
    totalNumSystems(ensemble.totalNumSystems),
    origin(ensemble.origin),
    invSteps(ensemble.invSteps),
    numSystems(ensemble.numSystems),
    type(ensemble.type),
    dimensionality(ensemble.dimensionality),
    centerPoint(ensemble.centerPoint)
{
    locks.resize(totalNumSystems);
    for (int i = 0; i < totalNumSystems; i++)
        omp_init_lock(&locks[i]);
    omp_init_lock(&logMutex);
    inParticlesLoop = false;
}


Ensemble& Ensemble::operator =(const Ensemble &ensemble)
{
    systems = ensemble.systems;
    pendingSystems = ensemble.pendingSystems;
    totalNumSystems = ensemble.totalNumSystems;
    origin = ensemble.origin;
    invSteps = ensemble.invSteps;
    numSystems = ensemble.numSystems;
    type = ensemble.type;
    inParticlesLoop = false;
    dimensionality = ensemble.dimensionality;
    centerPoint = ensemble.centerPoint;
    return *this;
}


Ensemble::~Ensemble()
{
    for (int i = 0; i < totalNumSystems; i++)
        omp_destroy_lock(&locks[i]);
    omp_destroy_lock(&logMutex);
}


void Ensemble::computeBoundaryIndexes()
{
    boundaryCellIndexes.clear();
    for (int i = 0; i < numSystems.x; ++i)
    for (int j = 0; j < numSystems.y; ++j)
    for (int k = 0; k < numSystems.z; ++k)
    {
        bool boundaryX = (i == 0 || i == numSystems.x - 1);
        bool boundaryY = (dimensionality > 1) && (j == 0 || j == numSystems.y - 1);
        bool boundaryZ = (dimensionality > 2) && (k == 0 || k == numSystems.z - 1);
        if (boundaryX || boundaryY || boundaryZ)
            boundaryCellIndexes.push_back(k + numSystems.z * (j + numSystems.y * i));
    }
}


void Ensemble::add(const Particle& newParticle)
{
    // workaround for modules adding new particles in 1D and 2D
    Particle particle = newParticle;
    for (int d = dimensionality; d < 3; d++)
        particle.coords[d] = centerPoint[d];

    Int3 idx = getSystemIdx(particle);
    if (!(idx >= Int3(0, 0, 0)) || !(idx < numSystems))
    {
        omp_set_lock(&logMutex);
        FP3 coords = particle.getPosition();
        FP3 momentum = particle.getMomentum();
        ///
        /*ComputationLog::getInstance().writeError("trying to add a particle outside of domain, coords = "
            + Utils::toString(coords) + ", momentum = " + Utils::toString(momentum) + ", type = " +
            Particle::types[particle.getType()].name());*/
        omp_unset_lock(&logMutex);
    }

    if (inParticlesLoop)
        pendingSystem(idx)->add(particle);
    else
    {
        ParticleSystem * s = system(idx);
        int linearIdx = idx.z + numSystems.z * (idx.y + numSystems.y * idx.x);
        omp_set_lock(&locks[linearIdx]);
        s->add(particle);
        omp_unset_lock(&locks[linearIdx]);
        if (s == currentSystem)
            currentParticle = currentSystem->currentParticle;
    }
}


void Ensemble::deleteCurrentParticle()
{
    currentSystem->deleteCurrentParticle();
    currentParticle = currentSystem->currentParticle;
}


void Ensemble::begin()
{
    addPendingParticles();
    inParticlesLoop = true;
    currentSystem = ptr(systems);
    while (currentSystem - ptr(systems) < totalNumSystems)
    {
        currentSystem->begin();
        if (!currentSystem->end())
        {
            currentParticle = currentSystem->currentParticle;
            return;
        }
        currentSystem++;
    }
}


void Ensemble::next()
{
    currentSystem->next();
    if (!currentSystem->end())
    {
        currentParticle = currentSystem->currentParticle;
    }
    else
    {
        currentSystem->begin();
        currentSystem++;
        while (currentSystem - ptr(systems) < totalNumSystems)
        {
            currentSystem->begin();
            if (!currentSystem->end())
            {
               currentParticle = currentSystem->currentParticle;
               return;
            }
            currentSystem++;
        }
    }
}


bool Ensemble::end()
{
    if (!currentSystem->end())
        return false;
    ParticleSystem * s = currentSystem + 1;
    while (s - ptr(systems) < totalNumSystems)
    {
        if (s->size())
            return false;
        s++;
    }
    inParticlesLoop = false;
    addPendingParticles();
    return true;
}


void Ensemble::addPendingParticles()
{
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < numSystems.x; i++)
    for (int j = 0; j < numSystems.y; j++)
    for (int k = 0; k < numSystems.z; k++)
    {
        ParticleSystem* src = pendingSystem(i, j, k);
        ParticleSystem* dst = system(i, j, k);
        for (int idx = 0; idx < src->size(); idx++)
            dst->add(src->particles[idx]);
        if (dst == currentSystem)
            currentParticle = dst->currentParticle;
        src->clear();
    }
}


void Ensemble::clear()
{
    for (int i = 0; i < totalNumSystems; ++i)
        systems[i].clear();
}


int Ensemble::getNumAllocatedParticles() const
{
    int result = 0;
    for (int i = 0; i < totalNumSystems; i++)
        result += (int)systems[i].particles.capacity() +
            (int)pendingSystems[i].particles.capacity();
    return result;
}


int Ensemble::getNumParticles() const
{
    int result = 0;
    for (int i = 0; i < totalNumSystems; i++)
        result += systems[i].size();
    return result;
}


void Ensemble::dumpParticles(Particle* particles, const Int3* minCellIdx,
    const Int3* maxCellIdx) const
{
    int particleIdx = 0;
    for (int i = minCellIdx->x; i < maxCellIdx->x; i++)
    for (int j = minCellIdx->y; j < maxCellIdx->y; j++)
    for (int k = minCellIdx->z; k < maxCellIdx->z; k++)
    {
        const ParticleSystem * s = system(i, j, k);
        for (int idx = 0; idx < s->size(); ++idx)
            particles[particleIdx++] = s->particles[idx];
    }
}


void Ensemble::dumpParticleInCellCount(int * particleInCellCount,
    const Int3 * minCellIdx, const Int3 * maxCellIdx) const
{
    int idx = 0;
    for (int i = minCellIdx->x; i < maxCellIdx->x; i++)
    for (int j = minCellIdx->y; j < maxCellIdx->y; j++)
    for (int k = minCellIdx->z; k < maxCellIdx->z; k++)
        particleInCellCount[idx++] = system(i, j, k)->size();
}


void Ensemble::dumpOutgoingParticles(std::vector<Particle>& particles)
{
    int numThreads = omp_get_max_threads();
    std::vector<int> particleCount(numThreads + 1, 0);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < boundaryCellIndexes.size(); i++)
        particleCount[omp_get_thread_num() + 1] += systems[boundaryCellIndexes[i]].size();

    for (int i = 1; i <= numThreads; i++)
        particleCount[i] += particleCount[i - 1];
    particles.resize(particleCount[numThreads]);

    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < boundaryCellIndexes.size(); idx++)
    {
        ParticleSystem& currentSystem = systems[boundaryCellIndexes[idx]];
        int& threadDstIdx = particleCount[omp_get_thread_num()];
        for (int i = 0; i < currentSystem.size(); i++)
            particles[threadDstIdx + i] = currentSystem.particles[i];
        threadDstIdx += currentSystem.size();
        currentSystem.clear();
    }
}


void Ensemble::shrinkToFit()
{
    #pragma omp parallel for schedule(static, 1)
    for (int i = 0; i < totalNumSystems; i++)
    {
        std::vector<Particle>(systems[i].particles).swap(systems[i].particles);
        systems[i].begin();
    }
}


void Ensemble::tryShrink()
{
    const double overuseThreshold = 2.5;
    double memoryOveruse = (double)getNumAllocatedParticles() / (double)getNumParticles();
    if (memoryOveruse < overuseThreshold)
        return;

    #pragma omp parallel for schedule(static, 1)
    for (int i = 0; i < totalNumSystems; i++)
    {
        std::vector<Particle>& particles = systems[i].particles;
        size_t maxAllowedSize = (size_t)(particles.size() * overuseThreshold);
        if (maxAllowedSize < particles.capacity())
        {
            // shrink-to-fit
            std::vector<Particle>(particles).swap(particles);
            systems[i].begin();
        }
        std::vector<Particle>& pendingParticles = pendingSystems[i].particles;
        size_t maxAllowedPendingSize = std::max(pendingParticles.size(), systems[i].particles.size());
        if (maxAllowedPendingSize < pendingParticles.capacity())
        {
            // shrink-to-fit
            std::vector<Particle>(pendingParticles).swap(pendingParticles);
            pendingSystems[i].begin();
        }
    }
}


std::string toUpper(const std::string& s)
{
    std::string result = s;
    for (int i = 0; i < result.length(); i++)
        result[i] = ::toupper(result[i]);
    return result;
}


int Ensemble::getTypeIdx(const std::string& name) const
{
    // First look for case-sensitive match
    for (int i = 0; i != type.size(); i++)
        if (type[i].name() == name)
            return i;

    // If it is not found, look for case-insensitive
    std::string nameUpper = toUpper(name);
    for (int i = 0; i != type.size(); i++)
    {
        std::string typeNameUpper = toUpper(type[i].name());
        if (typeNameUpper == nameUpper)
            return i;
    }

    ///
    ///ComputationLog::getInstance().writeError("requested particle type '" + name + "' is not defined");
    return -1;
}

} // namespace pica
