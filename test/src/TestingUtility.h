#ifndef TESTINGUTILITY_H
#define TESTINGUTILITY_H


#include "pica/Parameters.h"
#include "pica/Particle.h"
#include "pica/ParticleSystem.h"
#include "pica/Vectors.h"

#include "gtest/gtest.h"

#include <vector>


#define ASSERT_EQ_FP3(expected, actual) \
    ASSERT_EQ(expected.x, actual.x); \
    ASSERT_EQ(expected.y, actual.y); \
    ASSERT_EQ(expected.z, actual.z); \

#define ASSERT_EQ_INT3(expected, actual) \
    ASSERT_EQ(expected.x, actual.x); \
    ASSERT_EQ(expected.y, actual.y); \
    ASSERT_EQ(expected.z, actual.z); \

// Assert two FP3s are nearly equal:
// if expected value is not near zero, expect relative error is smaller than
// m_maxRelativeError, else expect absolute error is smaller than
// m_maxAbsoluteError.
#define ASSERT_NEAR_FP3(expected, actual) \
    if (expected.norm() > maxAbsoluteError) \
        ASSERT_LE(dist(expected, actual) / expected.norm(), maxRelativeError); \
    else \
        ASSERT_LE(dist(expected, actual), maxAbsoluteError); \

#define ASSERT_NEAR_FP(expected, actual) \
    if (fabs(expected) > maxAbsoluteError) \
        ASSERT_LE(fabs(expected - actual) / fabs(expected), maxRelativeError); \
    else \
        ASSERT_LE(fabs(expected - actual), maxAbsoluteError); \


class BaseFixture : public testing::Test {
protected:

    virtual void SetUp();
    virtual void TearDown();

    // Return whether two FP3s have coords differ by not larger than eps each.
    bool nearFP3(const pica::FP3 & a, const pica::FP3 & b, const pica::FP eps);
    // Return whether two particles have equal position, momentum and type index.
    bool eqParticle(const pica::Particle & a, const pica::Particle & b);

    // Get uniformly distributed in [a, b) pseudo-random number.
    pica::FP urand(pica::FP a, pica::FP b);
    int urandInt(int a, int b);
    // Get distributed in [a, b) pseudo-random vector.
    pica::FP3 urandFP3(pica::FP3 a, pica::FP3 b);
    pica::Int3 urandInt3(pica::Int3 a, pica::Int3 b);
    // Get _n_ random vectors between _minValue_ and _maxValue_.
    std::vector<pica::FP3> randomVectors(int n, pica::FP3 & minValue, pica::FP3 & maxValue);
    // Get point inside computational area.
    pica::FP3 internalPoint(const pica::Parameters& parameters);

    pica::FP maxAbsoluteError; // max absolute error that counts as "passed"
    pica::FP maxRelativeError; // max relative error that counts as "passed"
};


class BaseGridFixture : public BaseFixture {
protected:
    virtual void SetUp();
    pica::Parameters parameters;
};


class BaseParticleFixture : public BaseFixture {
protected:

    virtual void SetUp();
    pica::Particle randomParticle();

    pica::Parameters parameters;
};

bool eqParticles(const pica::Particle& a, const pica::Particle& b);
bool eqParticleSystems(pica::ParticleSystem& a, pica::ParticleSystem& b);


#endif
