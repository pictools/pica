#include "pica/ParticlePush.h"

#include "pica/currentDeposition/CurrentDeposition.h"
#include "pica/fieldInterpolation/TileInterpolator.h"
#include "pica/grid/Grid.h"
#include "pica/Ensemble.h"
#include "pica/threading/OpenMPHelper.h"
#include "pica/Parameters.h"
#include "pica/ParticleSystem.h"

#include <algorithm>


namespace pica {

void ParticlePush::init(const Parameters& _parameters, Ensemble& _ensemble, Grid& _grid)
{
    parameters = _parameters;
    ensemble = &_ensemble;
    grid = &_grid;
    cdt = constants::c * parameters.timeStep;
    coeff.resize(Particle::numTypes);
    for (int i = 0; i < Particle::numTypes; ++i)
        coeff[i] = Particle::types[i].charge * parameters.timeStep /
            ((FP)2 * Particle::types[i].mass * constants::c);
    centerPoint = 0.5 * (parameters.globalMin + parameters.globalMax);

    updateDims();
    
    // Choose Boris or Vay pusher
    useBoris = true;
    std::string type = "boris";
    for (int i = 0; i < type.length(); i++)
        type[i] = ::toupper(type[i]);
    if (type == "VAY")
        useBoris = false;
}


void ParticlePush::updateDims()
{
    Int3 numSystems = ensemble->numSystems;
    boundaryCellIndexes.clear();
    for (int i = 1; i < numSystems.x - 1; ++i)
    for (int j = 1; j < numSystems.y - 1; ++j)
    for (int k = 1; k < numSystems.z - 1; ++k)
    {
        if (isBoundaryCell(i, j, k))
            boundaryCellIndexes.push_back(Int3(i, j, k));
    }
}

void ParticlePush::runInner()
{
    if (useBoris)
        particlesLoopInner<BorisFunctor>();
    else
        particlesLoopInner<VayFunctor>();
}

void ParticlePush::runOuter()
{
    if (useBoris)
        particlesLoopOuter<BorisFunctor>();
    else
        particlesLoopOuter<VayFunctor>();
}

template <typename PusherType>
void ParticlePush::particlesLoopInner()
{
    switch (grid->getInterpolationType())
    {
        case Grid::Interpolation_CIC:
            particlesLoopInner<TileInterpolatorCIC, PusherType>();
            break;
        case Grid::Interpolation_TSC:
            if (parameters.globalGridSize.z > 1)
                particlesLoopInner<TileInterpolatorTSC, PusherType>();
            else
                if (parameters.globalGridSize.y > 1)
                    particlesLoopInner<XYTileInterpolatorTSC, PusherType>();
                else
                    particlesLoopInner<XTileInterpolatorTSC, PusherType>();
            break;
        default:
            particlesLoopInner<TileInterpolator<1>, PusherType>();
    }
}


template <typename TileInterpolatorType, typename PusherType>
void ParticlePush::particlesLoopInner()
{
    ensemble->inParticlesLoop = true;
    const Int3 numTiles = ensemble->numSystems;
    for (int ii = 0; ii < 3; ii++)
    for (int jj = 0; jj < 3; jj++)
    for (int kk = 0; kk < 3; kk++)
    {
        #pragma omp parallel for collapse(3) schedule(static, 1)
        for (int i = 1 + ii; i < numTiles.x - 1; i += 3)
        for (int j = 1 + jj; j < numTiles.y - 1; j += 3)
        for (int k = 1 + kk; k < numTiles.z - 1; k += 3)
        {
            ParticleSystem* system = ensemble->system(i, j, k);
            if (system->size() == 0)
                continue;

            if (!isBoundaryCell(i, j, k))
                particlesLoopInTile<TileInterpolatorType, PusherType>(i, j, k);

            handleMigration(i, j, k);
        }
    }

    ensemble->inParticlesLoop = false;
    ensemble->addPendingParticles();
}

template <typename PusherType>
void ParticlePush::particlesLoopOuter()
{
    switch (grid->getInterpolationType())
    {
        case Grid::Interpolation_CIC:
            particlesLoopOuter<TileInterpolatorCIC, PusherType>();
            break;
        case Grid::Interpolation_TSC:
            if (parameters.globalGridSize.z > 1)
                particlesLoopOuter<TileInterpolatorTSC, PusherType>();
            else
                if (parameters.globalGridSize.y > 1)
                    particlesLoopOuter<XYTileInterpolatorTSC, PusherType>();
                else
                    particlesLoopOuter<XTileInterpolatorTSC, PusherType>();
            break;
        default:
            particlesLoopOuter<TileInterpolator<1>, PusherType>();
    }
}

template <typename TileInterpolatorType, typename PusherType>
void ParticlePush::particlesLoopOuter()
{
    Int3 numSystems = ensemble->numSystems;
    ensemble->inParticlesLoop = true;

    #pragma omp parallel for
    for (int idx = 0; idx < (int)boundaryCellIndexes.size(); idx++)
    {
        const Int3& cell = boundaryCellIndexes[idx];
        particlesLoopInTile<TileInterpolatorType, PusherType>(cell.x, cell.y, cell.z);
    }
    ensemble->inParticlesLoop = false;

    #pragma omp parallel for
    for (int idx = 0; idx < (int)boundaryCellIndexes.size(); idx++)
    {
        const Int3& cell = boundaryCellIndexes[idx];
        ParticleSystem* system = ensemble->system(cell);
        for (int i = 0; i < system->size(); ++i)
        {
            const Particle& particle = system->particles[i];
            Int3 newCell = ensemble->getCellIndex(particle);
            if (newCell.x == 0 || newCell.x == numSystems.x - 1
                || newCell.y == 0 || newCell.y == numSystems.y - 1
                || newCell.z == 0 || newCell.z == numSystems.z - 1)
            {
                ensemble->add(particle);
                system->deleteParticle(i);
                --i;
            }
        }
    }

    ensemble->addPendingParticles();
}

template <typename TileInterpolatorType, typename PusherType>
void ParticlePush::particlesLoopInTile(int i, int j, int k)
{
    ParticleSystem* system = ensemble->system(i, j, k);
    TileInterpolatorType interpolator(i, j, k, grid);

    // Make interpolation and push for particles in chunks of fixed size.
    const int chunkSize = 32;
    Particle* particles = system->raw();
    const int numParticles = system->size();
    const int numChunks = (numParticles + chunkSize - 1) / chunkSize;
    for (int chunkIdx = 0; chunkIdx < numChunks; ++chunkIdx)
    {
        Particle* chunkParticles = particles + chunkIdx * chunkSize;
        int endIdx = std::min(chunkSize,
            system->size() - chunkIdx * chunkSize);

        FP3 interpolatedE[chunkSize], interpolatedB[chunkSize];
        FP eCoeff[chunkSize];
        #pragma novector
        for (int idx = 0; idx < endIdx; ++idx)
        {
            eCoeff[idx] = coeff[chunkParticles[idx].typeIndex];
            interpolator.getFields(chunkParticles[idx].coords,
                interpolatedE[idx], interpolatedB[idx]);
        }

        PusherType pusher(this);
        #pragma simd
        #pragma forceinline
        for (int idx = 0; idx < endIdx; ++idx)
            pusher(&chunkParticles[idx], interpolatedE[idx], interpolatedB[idx], eCoeff[idx]);
    }

    handleLowDimensionality(*system);
}

void ParticlePush::vay(Particle* particle, const FP3& e, const FP3& b, FP eCoeff)
{
    // The code below uses precomputed coefficient:
    // eCoeff = q * dt / (2 * m * c)
    FP3 eMomentum = e * eCoeff;
    FP3 tau = b * eCoeff;
    FP taunorm2 = tau.norm2();
    FP3 vbMomentum = cross(particle->p * particle->invGamma, tau);
    FP3 um = particle->p + eMomentum * (FP)2 + vbMomentum;

    FP u = dot(um, tau);
    FP sigma = (FP)1 + um.norm2() - taunorm2;
    FP gamma = sqrt((sigma + sqrt(sqr(sigma) +
        (FP)4 * (taunorm2 + u * u))) * (FP)0.5);
    particle->invGamma = (FP)1 / gamma;
    FP3 t = tau * particle->invGamma;
    FP s = (FP)1 / ((FP)1 + t.norm2());
    particle->p = s * (um + dot(um, t) * t + cross(um, t));
    particle->coords += particle->p * particle->invGamma * cdt;
}

void ParticlePush::boris(Particle* particle, const FP3& e, const FP3& b, FP eCoeff)
{
    // The code below uses precomputed coefficient:
    // eCoeff = q * dt / (2 * m * c)
    FP3 eMomentum = e * eCoeff;
    FP3 um = particle->p + eMomentum;
    FP3 t = b * eCoeff / sqrt((FP)1 + um.norm2());
    FP3 uprime = um + cross(um, t);
    FP3 s = t * (FP)2 / ((FP)1 + t.norm2());
    particle->p = eMomentum + um;
    particle->p += cross(uprime, s);
    particle->invGamma = (FP)1 / sqrt((FP)1 + particle->p.norm2());
    particle->coords += particle->p * particle->invGamma * cdt;
}

void ParticlePush::handleLowDimensionality(ParticleSystem& system)
{
    if (parameters.dimensionality < 3)
    {
        Particle* particles = system.raw();
        const int numParticles = system.size();
        // Fix coordinates along the ghost dimension to avoid migration along this dimension
        for (int d = parameters.dimensionality; d < 3; d++)
        {
            #pragma simd
            for (int idx = 0; idx < numParticles; idx++)
                particles[idx].coords[d] = centerPoint[d];
        }
    }
}


void ParticlePush::handleMigration(int i, int j, int k)
{
    ParticleSystem* system = ensemble->system(i, j, k);
    const Int3 tileIdx(i, j, k);
    for (int idx = 0; idx < system->size(); ++idx) {
        if (ensemble->isMigratingParticle(system->particles[idx], tileIdx))
        {
            ensemble->add(system->particles[idx]);
            system->deleteParticle(idx);
            --idx;
        }
    }
}


bool ParticlePush::isBoundaryCell(int i, int j, int k)
{
    Int3 numSystems = ensemble->numSystems;
    int dimensionality = parameters.dimensionality;
    bool boundaryX = (i == 1 || i == numSystems.x - 2);
    bool boundaryY = (j == 1 || j == numSystems.y - 2) && (dimensionality >= 2);
    bool boundaryZ = (k == 1 || k == numSystems.z - 2) && (dimensionality == 3);
    return boundaryX || boundaryY || boundaryZ;
}

} // namespace pica
