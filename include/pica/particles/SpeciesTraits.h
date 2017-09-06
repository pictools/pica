#ifndef PICA_SPECIES_TRAITS_H
#define PICA_SPECIES_TRAITS_H


#include "pica/particles/ParticleTraits.h"


namespace pica {


// This class defines types of species data
// By default it relies on nested types in a species class
// In case a species class does not have these, SpeciesTraits must be specialized
template<typename SpeciesType>
struct SpeciesTraits {
public:
    typedef SpeciesType::ChargeType ChargeType;
    typedef SpeciesType::MassType MassType;
};

// Specialization for the common case, use types from ParticleTraits
template<typename ParticleType>
struct SpeciesTraits<typename ParticleType::Species> {
public:
    typedef ParticleTraits<ParticleType>::ChargeType ChargeType;
    typedef ParticleTraits<ParticleType>::MassType MassType;
};


} // namespace pica


#endif

