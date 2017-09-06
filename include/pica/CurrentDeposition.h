#ifndef PICA_CURRENTDEPOSITION_H
#define PICA_CURRENTDEPOSITION_H


#include "pica/Parameters.h"
#include "pica/utility/Utility.h"
#include "pica/math/Vectors.h"

#include <memory>
#include <vector>


namespace pica {

class Ensemble;
class Grid;

class CurrentDeposition {
public:

    void init(const Parameters& parameters, Grid& grid);
    void run(const Ensemble& ensemble, Grid& grid);

    Grid* grid; ///
    Parameters parameters;
    std::vector<FP> currentDensityCoeff; // precomputed coefficients
    std::vector<FP3> esirkepovCurrentDensityCoeff; // precomputed coefficients
    std::vector<FP> chargeCoeff;

private:

    template <typename TileDepositorType>
    void addCurrents(const Ensemble& ensemble, Grid& grid);

    template <typename TileDepositorType>
    void addCurrentsFromCell(int i, int j, int k, const Ensemble& ensemble);

};

} // namespace pica


#endif
