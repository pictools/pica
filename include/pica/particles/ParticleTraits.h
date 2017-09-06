#ifndef PICA_PARTICLE_TRAITS_H
#define PICA_PARTICLE_TRAITS_H


namespace pica {


// This class defines types of particle data
// By default it relies on nested types in a particle class
// In case a particle class does not have these, ParticleTraits must be specialized
template<typename ParticleType>
struct ParticleTraits {
public:
    typedef ParticleType::PositionType PositionType;
    typedef ParticleType::MomemtumType MomemtumType;
    typedef ParticleType::GammaType GammaType;
    typedef ParticleType::ChargeType ChargeType;
    typedef ParticleType::MassType MassType;
    typedef ParticleType::FactorType FactorType;
    typedef ParticleType::SpeciesType SpeciesType;
};


} // namespace pica


#endif

