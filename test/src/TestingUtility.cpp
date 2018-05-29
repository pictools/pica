#include "TestingUtility.h"

#include "pica/math/Vectors.h"
#include "pica/particles/Particle.h"
#include "pica/utility/Utility.h"

#include <cmath>
#include <memory>

using namespace pica;

namespace pica {

    namespace ParticleTypes {
        std::vector<ParticleType> typesVector;
        const ParticleType* types = nullptr;
        int numTypes = 0;
    } // namespace ParticleTypes

} // namespace pica

void BaseFixture::SetUp() {
    srand(1);
    maxAbsoluteError = (FP)1e-4;
    maxRelativeError = (FP)1e-4;
}

void BaseFixture::TearDown() {
}

// Return whether two FP3s have coords differ by not larger than eps each.
bool BaseFixture::nearFP3(const FP3 & a, const FP3 & b, const FP eps) {
    return (fabs(a.x - b.x) <= eps) && (fabs(a.y - b.y) <= eps) &&
           (fabs(a.z - b.z) <= eps);
}

// Get uniformly distributed in [a, b) pseudo-random number.
FP BaseFixture::urand(FP a, FP b) const {
    return a + (b - a) * ((FP)rand()) / RAND_MAX;
}

int BaseFixture::urandInt(int a, int b) {
    return a + rand() % (b - a + 1);
}

// Get distributed in [a, b) pseudo-random vector.
FP3 BaseFixture::urandFP3(FP3 a, FP3 b) {
    FP3 result;
    result.x = urand(a.x, b.x);
    result.y = urand(a.y, b.y);
    result.z = urand(a.z, b.z);
    return result;
}

Int3 BaseFixture::urandInt3(Int3 a, Int3 b) {
    Int3 result;
    result.x = urandInt(a.x, b.x);
    result.y = urandInt(a.y, b.y);
    result.z = urandInt(a.z, b.z);
    return result;
}

// Get _n_ random vectors between _minValue_ and _maxValue_.
std::vector<FP3> BaseFixture::randomVectors(int n, FP3 & minValue, FP3 & maxValue) {
    std::vector<FP3> result(n);
    for (int i = 0; i < n; ++i)
        result[i] = urandFP3(minValue, maxValue);
    return result;
}

// Get point inside computational area.
FP3 BaseFixture::internalPoint(const Parameters& parameters) {
    return urandFP3(parameters.globalMin, parameters.globalMax);
}


void BaseGridFixture::SetUp() {
    BaseFixture::SetUp();
    maxAbsoluteError = (FP)1e-4;
    maxRelativeError = (FP)0.05;

    parameters.localMin = parameters.globalMin = FP3(1.30986, -15.92346, -0.45621);
    parameters.localMax = parameters.globalMax = FP3(2.969704, -13.01246, 2.2919123);
    parameters.localGridSize = parameters.globalGridSize = parameters.gridSize = Int3(36, 42, 39);
    parameters.localMinIndex = Int3(0, 0, 0);
    parameters.localMaxIndex = parameters.localGridSize - Int3(1, 1, 1);
    parameters.step = (parameters.globalMax - parameters.globalMin) / FP3(parameters.globalGridSize);
    parameters.numIterations = 1;
    parameters.timeStep = 1.463e-12;
}
