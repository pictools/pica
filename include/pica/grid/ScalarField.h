#ifndef PICA_SCALARFIELD_H
#define PICA_SCALARFIELD_H


#include "pica/math/Vectors.h"

#include <vector>

namespace pica {

/* Class for storing 3d scalar field on a regular grid.
Provides index-wise access, interpolation and deposition.*/
class ScalarField
{
public:

    ScalarField(const Int3& size);
    ScalarField(const ScalarField& field);
    ScalarField& operator =(const ScalarField& field);

    /* Read-only access by scalar indexes */
    FP operator()(int i, int j, int k) const
    { return raw[k + (j + i * size.y) * size.z]; }

    /* Read-write access by scalar indexes */
    FP& operator()(int i, int j, int k)
    { return raw[k + (j + i * size.y) * size.z]; }

    /* Read-only access by vector index */
    FP operator()(const Int3& index) const
    { return raw[index.z + (index.y + index.x * size.y) * size.z]; }

    /* Read-write access by vector index */
    FP& operator()(const Int3& index)
    { return raw[index.z + (index.y + index.x * size.y) * size.z]; }

    /* Interpolation: with given base index and coefficients */
    FP interpolateCIC(const Int3& baseIdx, const FP3& coeffs) const;
    FP interpolateTSC(const Int3& baseIdx, const FP3& coeffs) const;
    FP interpolateSecondOrder(const Int3& baseIdx, const FP3& coeffs) const;
    FP interpolateFourthOrder(const Int3& baseIdx, const FP3& coeffs) const;
    FP interpolatePCS(const Int3& baseIdx, const FP3& coeffs) const;

    /* Deposition with given base index and coefficients */
    void depositCIC(FP value, const Int3& baseIdx, const FP3& coeffs);
    void depositTSC(FP value, const Int3& baseIdx, const FP3& coeffs);
    void depositPCS(FP value, const Int3& baseIdx, const FP3& coeffs);

    /* Make all elements zero */
    void zeroize();

private:

    FP interpolateThreePoints(const Int3& baseIdx, FP c[3][3]) const;

    std::vector<FP> elements; // storage
    FP* raw; // raw pointer to elements vector
    Int3 size; // size of each dimension
    Int3 dimensionCoeffInt; // 0 for fake dimensions, 1 otherwise
    FP3 dimensionCoeffFP; // 0 for fake dimensions, 1 otherwise

};

} // namespace pica


#endif
