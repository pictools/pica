#ifndef PICA_BENCHMARK_ENSEMBLEHELPER_H
#define PICA_BENCHMARK_ENSEMBLEHELPER_H


#include "pica/math/Vectors.h"
#include "pica/particles/Ensemble.h"
#include "pica/particles/EnsembleSupercells.h"

#include <memory>


namespace internal {

template<class Ensemble>
struct EnsembleFactory {
    typedef typename Ensemble::PositionType PositionType;

    static std::auto_ptr<Ensemble> create(PositionType minPosition, PositionType maxPosition,
        pica::Int3, pica::Int3)
    {
        return std::auto_ptr<Ensemble>(new Ensemble(minPosition, maxPosition));
    }
};

template<class ParticleArray>
struct EnsembleFactory<pica::EnsembleSupercells<ParticleArray>> {
    typedef pica::EnsembleSupercells<ParticleArray> Ensemble;
    typedef typename Ensemble::PositionType PositionType;

    static std::auto_ptr<Ensemble> create(PositionType minPosition, PositionType maxPosition,
        pica::Int3 numCells, pica::Int3 numCellsPerSupercell)
    {
        return std::auto_ptr<Ensemble>(new Ensemble(minPosition, maxPosition, numCells, numCellsPerSupercell));
    }
};


}


template<class Ensemble>
std::auto_ptr<Ensemble> createEnsemble(typename Ensemble::PositionType minPosition,
    typename Ensemble::PositionType maxPosition, pica::Int3 numCells, pica::Int3 numCellsPerSupercell = pica::Int3())
{
    return ::internal::EnsembleFactory<Ensemble>::create(minPosition, maxPosition, numCells, numCellsPerSupercell);
}


#endif
