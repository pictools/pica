#ifndef PICA_PARTICLE_TRAITS_H
#define PICA_PARTICLE_TRAITS_H


namespace pica {


// This class defines types of particle data
// By default it relies on nested types in a particle class
// In case a particle class does not have these, ParticleTraits must be specialized
template<typename ParticleType>
struct ParticleTraits {
public:
    typedef typename ParticleType::PositionType PositionType;
    typedef typename ParticleType::MomentumType MomentumType;
    typedef typename ParticleType::GammaType GammaType;
    typedef typename ParticleType::FactorType FactorType;
    typedef typename ParticleType::TypeIndexType TypeIndexType;
};


} // namespace pica


#endif

