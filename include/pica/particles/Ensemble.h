#ifndef PICA_ENSEMBLE_H
#define PICA_ENSEMBLE_H


#include "pica/math/Dimension.h"
#include "pica/particles/EnsembleOrdered.h"
#include "pica/particles/EnsembleSupercells.h"
#include "pica/particles/EnsembleUnordered.h"
#include "pica/particles/Particle.h"
#include "pica/particles/ParticleTraits.h"
#include "pica/threading/OpenMPHelper.h"

#include <algorithm>
#include <map>
#include <string>
#include <vector>


namespace pica {


enum EnsembleRepresentation { EnsembleRepresentation_Unordered, EnsembleRepresentation_Ordered,
    EnsembleRepresentation_Supercells };

inline std::string toString(EnsembleRepresentation ensembleRepresentation)
{
    std::map<EnsembleRepresentation, std::string> names;
    names[EnsembleRepresentation_Unordered] = "unordered";
    names[EnsembleRepresentation_Ordered] = "ordered";
    names[EnsembleRepresentation_Supercells] = "supercells";
    return names[ensembleRepresentation];
}

// Traits class to provide a Type corresponding to array of particles
// according to the given representation
template<class ParticleArray, EnsembleRepresentation storage>
struct Ensemble {
};

template<class ParticleArray>
struct Ensemble<ParticleArray, EnsembleRepresentation_Unordered> {
    typedef EnsembleUnordered<ParticleArray> Type;
};

template<class ParticleArray>
struct Ensemble<ParticleArray, EnsembleRepresentation_Ordered> {
    typedef EnsembleOrdered<ParticleArray> Type;
};

template<class ParticleArray>
struct Ensemble<ParticleArray, EnsembleRepresentation_Supercells> {
    typedef EnsembleSupercells<ParticleArray> Type;
};


} // namespace pica


#endif
