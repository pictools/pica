#include "pica/grid/ScalarField.h"

#include "pica/Formfactor.h"
#include "pica/threading/OpenMPHelper.h"
#include "pica/utility/Utility.h"


namespace pica {

ScalarField::ScalarField(const Int3& _size)
{
    size = _size;
    elements.resize(size.volume());
    raw = ptr(elements);
    for (int d = 0; d < 3; d++) {
        dimensionCoeffInt[d] = (size[d] > 1) ? 1 : 0;
        dimensionCoeffFP[d] = (FP)dimensionCoeffInt[d];
    }
    zeroize();
}


ScalarField::ScalarField(const ScalarField& field)
{
    size = field.size;
    elements = field.elements;
    raw = ptr(elements);
    dimensionCoeffInt = field.dimensionCoeffInt;
    dimensionCoeffFP = field.dimensionCoeffFP;
}


ScalarField& ScalarField::operator =(const ScalarField& field)
{
    size = field.size;
    elements = field.elements;
    raw = ptr(elements);
    dimensionCoeffInt = field.dimensionCoeffInt;
    dimensionCoeffFP = field.dimensionCoeffFP;
    return *this;
}


FP ScalarField::interpolateCIC(const Int3& baseIdx, const FP3& coeffs) const
{
    FP3 c = coeffs * dimensionCoeffFP;
    FP3 invC = FP3(1, 1, 1) - c;
    Int3 base = baseIdx * dimensionCoeffInt;
    Int3 next = base + dimensionCoeffInt;
    return invC.x * (invC.y * (invC.z * (*this)(base.x, base.y, base.z) + c.z * (*this)(base.x, base.y, next.z)) +
                    c.y * (invC.z * (*this)(base.x, next.y, base.z) + c.z * (*this)(base.x, next.y, next.z))) +
        c.x * (invC.y * (invC.z * (*this)(next.x, base.y, base.z) + c.z * (*this)(next.x, base.y, next.z)) +
        c.y * (invC.z * (*this)(next.x, next.y, base.z) + c.z * (*this)(next.x, next.y, next.z)));
}


FP ScalarField::interpolateTSC(const Int3& baseIdx, const FP3& coeffs) const
{
    FP c[3][3];
    for (int i = 0; i < 3; i++)
        c[0][i] = formfactorTSC(FP(i - 1) - coeffs.x);
    for (int j = 0; j < 3; j++)
        c[1][j] = formfactorTSC(FP(j - 1) - coeffs.y);
    for (int k = 0; k < 3; k++)
        c[2][k] = formfactorTSC(FP(k - 1) - coeffs.z);
    return interpolateThreePoints(baseIdx, c);
}


FP ScalarField::interpolateSecondOrder(const Int3& baseIdx, const FP3& coeffs) const
{
    FP c[3][3];
    c[0][0] = (FP)0.5 * (coeffs.x * (coeffs.x - (FP)1));
    c[0][1] = (FP)1 - coeffs.x * coeffs.x;
    c[0][2] = (FP)0.5 * (coeffs.x * (coeffs.x + (FP)1));
    c[1][0] = (FP)0.5 * (coeffs.y * (coeffs.y - (FP)1));
    c[1][1] = (FP)1 - coeffs.y * coeffs.y;
    c[1][2] = (FP)0.5 * (coeffs.y * (coeffs.y + (FP)1));
    c[2][0] = (FP)0.5 * (coeffs.z * (coeffs.z - (FP)1));
    c[2][1] = (FP)1 - coeffs.z * coeffs.z;
    c[2][2] = (FP)0.5 * (coeffs.z * (coeffs.z + (FP)1));
    return interpolateThreePoints(baseIdx, c);
}


FP ScalarField::interpolateThreePoints(const Int3& baseIdx, FP c[3][3]) const
{
    for (int d = 0; d < 3; d++)
        if (!dimensionCoeffInt[d]) {
            c[d][0] = 0;
            c[d][1] = 1.0;
            c[d][2] = 0;
        }
    FP result = 0;
    Int3 minIndex = Int3(-1, -1, -1) * dimensionCoeffInt;
    Int3 maxIndex = Int3(1, 1, 1) * dimensionCoeffInt;
    Int3 base = baseIdx * dimensionCoeffInt;
    for (int ii = minIndex.x; ii <= maxIndex.x; ii++)
    for (int jj = minIndex.y; jj <= maxIndex.y; jj++)
    for (int kk = minIndex.z; kk <= maxIndex.z; kk++)
        result += c[0][ii + 1] * c[1][jj + 1] * c[2][kk + 1] * (*this)(base.x + ii, base.y + jj, base.z + kk);
    return result;
}


FP ScalarField::interpolateFourthOrder(const Int3& baseIdx, const FP3& coeffs) const
{
    Int3 base = baseIdx * dimensionCoeffInt;
    const Int3 minAllowedIdx = Int3(2, 2, 2) * dimensionCoeffInt;
    const Int3 maxAllowedIdx = (size - Int3(4, 4, 4)) * dimensionCoeffInt;
    if ((base >= minAllowedIdx) && (base <= maxAllowedIdx) == false)
        return interpolateTSC(baseIdx, coeffs);
    FP c[3][5];
    formfactorFourthOrder(coeffs.x, c[0]);
    formfactorFourthOrder(coeffs.y, c[1]);
    formfactorFourthOrder(coeffs.z, c[2]);
    for (int d = 0; d < 3; d++)
        if (!dimensionCoeffInt[d]) {
            c[d][0] = 0;
            c[d][1] = 0;
            c[d][2] = 1.0;
            c[d][3] = 0;
            c[d][4] = 0;
        }
    FP result = 0;
    Int3 minIndex = Int3(-2, -2, -2) * dimensionCoeffInt;
    Int3 maxIndex = Int3(2, 2, 2) * dimensionCoeffInt;
    for (int ii = minIndex.x; ii <= maxIndex.x; ii++)
    for (int jj = minIndex.y; jj <= maxIndex.y; jj++)
    for (int kk = minIndex.z; kk <= maxIndex.z; kk++)
        result += c[0][ii + 2] * c[1][jj + 2] * c[2][kk + 2] * (*this)(base.x + ii, base.y + jj, base.z + kk);
    return result;
}


FP ScalarField::interpolatePCS(const Int3& baseIdx, const FP3& coeffs) const
{
    FP c[3][4];
    for (int i = 0; i < 4; i++)
        c[0][i] = formfactorPCS(FP(i - 1) - coeffs.x);
    for (int j = 0; j < 4; j++)
        c[1][j] = formfactorPCS(FP(j - 1) - coeffs.y);
    for (int k = 0; k < 4; k++)
        c[2][k] = formfactorPCS(FP(k - 1) - coeffs.z);
    for (int d = 0; d < 3; d++)
        if (!dimensionCoeffInt[d]) {
            c[d][0] = 1.0;
            c[d][1] = 0;
            c[d][2] = 0;
            c[d][3] = 0;
        }
    FP result = 0;
    Int3 base = (baseIdx - Int3(1, 1, 1)) * dimensionCoeffInt;
    Int3 minIndex = Int3(0, 0, 0);
    Int3 maxIndex = Int3(3, 3, 3) * dimensionCoeffInt;
    for (int ii = minIndex.x; ii <= maxIndex.x; ii++)
    for (int jj = minIndex.y; jj <= maxIndex.y; jj++)
    for (int kk = minIndex.z; kk <= maxIndex.z; kk++)
        result += c[0][ii] * c[1][jj] * c[2][kk] * (*this)(base.x + ii, base.y + jj, base.z + kk);
    return result;
}


void ScalarField::depositCIC(FP value, const Int3& baseIdx, const FP3& coeffs)
{
    FP3 c = coeffs;
    FP3 invC = FP3(1, 1, 1) - c;
    Int3 base = baseIdx * dimensionCoeffInt;
    Int3 next = base + dimensionCoeffInt;
    (*this)(base.x, base.y, base.z) += invC.x * invC.y * invC.z * value;
    (*this)(base.x, base.y, next.z) += invC.x * invC.y * c.z * value;
    (*this)(base.x, next.y, base.z) += invC.x * c.y * invC.z * value;
    (*this)(base.x, next.y, next.z) += invC.x * c.y * c.z * value;
    (*this)(next.x, base.y, base.z) += c.x * invC.y * invC.z * value;
    (*this)(next.x, base.y, next.z) += c.x * invC.y * c.z * value;
    (*this)(next.x, next.y, base.z) += c.x * c.y * invC.z * value;
    (*this)(next.x, next.y, next.z) += c.x * c.y * c.z * value;
}


void ScalarField::depositTSC(FP value, const Int3& baseIdx, const FP3& coeffs)
{
    FP c[3][3];
    for (int i = 0; i < 3; i++)
        c[0][i] = formfactorTSC(FP(i - 1) - coeffs.x);
    for (int j = 0; j < 3; j++)
        c[1][j] = formfactorTSC(FP(j - 1) - coeffs.y);
    for (int k = 0; k < 3; k++)
        c[2][k] = formfactorTSC(FP(k - 1) - coeffs.z);
    Int3 minIndex = Int3(-1, -1, -1) * dimensionCoeffInt;
    Int3 maxIndex = Int3(1, 1, 1) * dimensionCoeffInt;
    for (int d = 0; d < 3; d++)
        if (!dimensionCoeffInt[d]) {
            c[d][0] = 0;
            c[d][1] = 1.0;
            c[d][2] = 0;
        }
    FP result = 0;
    Int3 base = baseIdx * dimensionCoeffInt;
    for (int ii = minIndex.x; ii <= maxIndex.x; ii++)
    for (int jj = minIndex.y; jj <= maxIndex.y; jj++)
    for (int kk = minIndex.z; kk <= maxIndex.z; kk++)
        (*this)(base.x + ii, base.y + jj, base.z + kk) += c[0][ii + 1] * c[1][jj + 1] * c[2][kk + 1] * value;
}


void ScalarField::depositPCS(FP value, const Int3& baseIdx, const FP3& coeffs)
{
    FP c[3][4];
    for (int i = 0; i < 4; i++)
        c[0][i] = formfactorPCS(FP(i - 1) - coeffs.x);
    for (int j = 0; j < 4; j++)
        c[1][j] = formfactorPCS(FP(j - 1) - coeffs.y);
    for (int k = 0; k < 4; k++)
        c[2][k] = formfactorPCS(FP(k - 1) - coeffs.z);
    for (int d = 0; d < 3; d++)
        if (!dimensionCoeffInt[d]) {
            c[d][0] = 1.0;
            c[d][1] = 0;
            c[d][2] = 0;
            c[d][3] = 0;
        }
    FP result = 0;
    Int3 base = (baseIdx - Int3(1, 1, 1)) * dimensionCoeffInt;
    Int3 minIndex = Int3(0, 0, 0);
    Int3 maxIndex = Int3(3, 3, 3) * dimensionCoeffInt;
    for (int ii = minIndex.x; ii <= maxIndex.x; ii++)
    for (int jj = minIndex.y; jj <= maxIndex.y; jj++)
    for (int kk = minIndex.z; kk <= maxIndex.z; kk++)
        (*this)(base.x + ii, base.y + jj, base.z + kk) += c[0][ii] * c[1][jj] * c[2][kk] * value;
}


void ScalarField::zeroize()
{
    #pragma omp parallel for
    for (int idx = 0; idx < (int)elements.size(); idx++)
        elements[idx] = 0;
}

} // namespace pica
