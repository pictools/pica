#ifndef PICA_PARTICLEPUSH_H
#define PICA_PARTICLEPUSH_H


#include "pica/Grid.h"
#include "pica/Ensemble.h"
#include "pica/Parameters.h"
#include "pica/particle/Particle.h"
#include "pica/ParticleSystem.h"
#include "pica/math/Vectors.h"
#include "pica/utility/Utility.h"


namespace pica {

// The particle push stage of PIC. Moves all particles and leaves ensemble in
// correct state.
class ParticlePush {

    friend class BorisFunctor;
    friend class VayFunctor;

public:

    void init(const Parameters& _parameters, Ensemble& ensemble, Grid& grid);
    void runInner();
    void runOuter();

    void updateDims();

private:

    Parameters parameters;
    Ensemble* ensemble;
    Grid* grid;
    std::vector<FP> coeff; // precomputed coefficients for Boris push
    FP cdt;
    FP3 centerPoint; // for 1D and 2D simulations
    std::vector<Int3> boundaryCellIndexes;
    bool useBoris;

    template <typename PusherType>
    void particlesLoopInner();
    template <typename TileInterpolatorType, typename PusherType>
    void particlesLoopInner();

    template <typename PusherType>
    void particlesLoopOuter();
    template <typename TileInterpolatorType, typename PusherType>
    void particlesLoopOuter();

    template <typename TileInterpolatorType, typename PusherType>
    void particlesLoopInTile(int i, int j, int k);
    void boris(Particle* particle, const FP3& e, const FP3& b, FP eCoeff);
    void vay(Particle* particle, const FP3& e, const FP3& b, FP eCoeff);
    void handleLowDimensionality(ParticleSystem& system);
    void handleMigration(int i, int j, int k);
    bool isBoundaryCell(int i, int j, int k);
};


class BorisFunctor {
public:
    BorisFunctor(ParticlePush* _particlePush): particlePush(_particlePush) {} 
    inline void operator() (Particle* particle, const FP3& e, const FP3& b, FP eCoeff) {
        particlePush->boris(particle, e, b, eCoeff);
    }
private:
    ParticlePush* particlePush;
};

class VayFunctor {
public:
    VayFunctor(ParticlePush* _particlePush): particlePush(_particlePush) {} 
    inline void operator() (Particle* particle, const FP3& e, const FP3& b, FP eCoeff) {
        particlePush->vay(particle, e, b, eCoeff);
    }
private:
    ParticlePush* particlePush;
};

} // namespace pica


#endif
