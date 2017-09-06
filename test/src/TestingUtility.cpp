#include "TestingUtility.h"

#include "pica/Particle.h"
#include "pica/Utility.h"
#include "pica/math/Vectors.h"

#include <cmath>
#include <memory>

using namespace pica;


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

// Return whether two particles have equal position, momentum and type index.
bool BaseFixture::eqParticle(const Particle & a, const Particle & b) {
    // since coords are transformed inside solver we don't require
    // eqFP3 on coords but only nearFP3
    return nearFP3(a.getPosition(), b.getPosition(), 1e-5) &&
           (a.getMomentum() == b.getMomentum()) && (a.getType() == b.getType());
}

// Get uniformly distributed in [a, b) pseudo-random number.
FP BaseFixture::urand(FP a, FP b) {
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
    parameters.step = (parameters.globalMax - parameters.globalMin) / parameters.globalGridSize;
    parameters.numIterations = 1;
    parameters.timeStep = 1.463e-12;
}


void BaseParticleFixture::SetUp() {
    BaseFixture::SetUp();
    maxAbsoluteError = (FP)1e-12;
    maxRelativeError = (FP)0.005; // 0.5%

    parameters.localMin = parameters.globalMin = FP3(-6.8913, -3.2154, 2.0914);
    parameters.localMax = parameters.globalMax = FP3(9.6214, -0.9123, 8.2135);
    parameters.localGridSize = parameters.globalGridSize = parameters.gridSize = Int3(4, 7, 6);
    parameters.localMinIndex = Int3(0, 0, 0);
    parameters.localMaxIndex = parameters.localGridSize - Int3(1, 1, 1);
    parameters.step = (parameters.globalMax - parameters.globalMin) / parameters.globalGridSize;
    parameters.timeStep = 1e-12;
    parameters.numIterations = 1;

    Particle::typesVector.resize(2);
    Particle::typesVector[0].mass = constants::electronMass;
    Particle::typesVector[0].charge = constants::electronCharge;
    Particle::typesVector[0].name() = "Electron";
    Particle::typesVector[1].mass = constants::protonMass;
    Particle::typesVector[1].charge = -constants::electronCharge;
    Particle::typesVector[1].name() = "Proton";
    Particle::types = ptr(Particle::typesVector);
    Particle::numTypes = (int)Particle::typesVector.size();
}

Particle BaseParticleFixture::randomParticle() {
    FP minCoeff = (FP)-10;
    FP maxCoeff = (FP)10;
    FP minMomentum = minCoeff;
    FP maxMomentum = maxCoeff;
    FP3 momentum(urand(minMomentum, maxMomentum),
        urand(minMomentum, maxMomentum),
        urand(minMomentum, maxMomentum));
    float factor = (float)urand(1e-5, 1e5);
    return Particle(internalPoint(parameters), momentum, rand() % Particle::numTypes, factor);
}

bool eqParticles(const Particle& a, const Particle& b) {
    return (a.getPosition() == b.getPosition()) &&
        (a.getMomentum() == b.getMomentum()) &&
        (a.getFactor() == b.getFactor()) &&
        (a.getType() == b.getType());
}

bool eqParticleSystems(ParticleSystem& a, ParticleSystem& b) {
    if (a.size() != b.size())
        return false;
    a.begin();
    b.begin();
    for (; !a.end(); ) {
        if (!eqParticles(*a.currentParticle, *b.currentParticle))
            return false;
        a.next();
        b.next();
    }
    return true;
}
