#ifndef PICA_SPECIES_INTERFACE_H
#define PICA_SPECIES_INTERFACE_H

// This file serves as a documentation and should not be included as it will fail to build

namespace pica {


// This class defines an interface for particle classes to follow
class SpeciesInterface {
public:

    // Either provide the following types or specialization of SpeciesTraits
    typedef double ChargeType;
    typedef double MassType;

    MassType getMass() const;
    ChargeType getCharge() const;
    std::string getName() const;

};


} // namespace pica


#endif

