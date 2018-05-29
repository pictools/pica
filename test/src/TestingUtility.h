#ifndef TESTINGUTILITY_H
#define TESTINGUTILITY_H


#include "pica/Parameters.h"
#include "pica/math/Constants.h"
#include "pica/math/Vectors.h"
#include "pica/particles/Particle.h"
#include "pica/particles/ParticleTraits.h"

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

#define ASSERT_EQ_VECTOR(expected, actual, dimension) \
    for (int d = 0; d < dimension; d++) \
        ASSERT_EQ(expected[d], actual[d]); \

// Assert two FP3s are nearly equal:
// if expected value is not near zero, expect relative error is smaller than
// m_maxRelativeError, else expect absolute error is smaller than
// m_maxAbsoluteError.
#define ASSERT_NEAR_VECTOR(expected, actual) \
    if (expected.norm() > this->maxAbsoluteError) \
        ASSERT_LE(dist(expected, actual) / expected.norm(), this->maxRelativeError); \
    else \
        ASSERT_LE(dist(expected, actual), this->maxAbsoluteError);

#define ASSERT_NEAR_FP3(expected, actual) \
    if (expected.norm() > maxAbsoluteError) \
        ASSERT_LE(dist(expected, actual) / expected.norm(), maxRelativeError); \
    else \
        ASSERT_LE(dist(expected, actual), maxAbsoluteError); \

#define ASSERT_NEAR_FP(expected, actual) \
    if (fabs(expected) > this->maxAbsoluteError) \
        ASSERT_LE(fabs(expected - actual) / fabs(expected), this->maxRelativeError); \
    else \
        ASSERT_LE(fabs(expected - actual), this->maxAbsoluteError); \


class BaseFixture : public testing::Test {
protected:

    virtual void SetUp();
    virtual void TearDown();

    // Return whether two FP3s have coords differ by not larger than eps each.
    bool nearFP3(const pica::FP3 & a, const pica::FP3 & b, const pica::FP eps);

    // Get uniformly distributed in [a, b) pseudo-random number.
    pica::FP urand(pica::FP a, pica::FP b) const;
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


template<class ParticleType>
class BaseParticleFixture : public BaseFixture {
public:
    typedef ParticleType Particle;
    typedef typename pica::ParticleTraits<Particle>::PositionType PositionType;
    typedef typename pica::ParticleTraits<Particle>::MomentumType MomentumType;
    typedef typename pica::ParticleTraits<Particle>::GammaType GammaType;
    typedef typename pica::ParticleTraits<Particle>::FactorType FactorType;
    typedef typename pica::ScalarType<MomentumType>::Type Real;
    static const int dimension = pica::VectorDimensionHelper<PositionType>::dimension;
    static const int momentumDimension = pica::VectorDimensionHelper<MomentumType>::dimension;

    virtual void SetUp()
    {
        BaseFixture::SetUp();

        pica::ParticleTypes::numTypes = 1;
        pica::ParticleTypes::typesVector.resize(1);
        pica::ParticleTypes::types = &pica::ParticleTypes::typesVector[0];
        pica::ParticleTypes::typesVector[0].mass = pica::Constants<pica::MassType>::electronMass();
        pica::ParticleTypes::typesVector[0].charge = pica::Constants<pica::ChargeType>::electronCharge();
    }

    // Helper function to unify initialization of positions for 1d, 2d and 3d
    // In 1d y, z are ignored, in 2d z is ignored
    PositionType getPosition(Real x, Real y, Real z) const
    {
        Real positionArray[] = { x, y, z };
        PositionType position;
        for (int d = 0; d < dimension; d++)
            position[d] = positionArray[d];
        return position;
    }

    Particle randomParticle()
    {
        Real minPosition = -10;
        Real maxPosition = 10;
        return randomParticle(getPosition(minPosition, minPosition, minPosition),
            getPosition(maxPosition, maxPosition, maxPosition));
    }

    Particle randomParticle(PositionType minPosition, PositionType maxPosition)
    {
        PositionType position;
        for (int d = 0; d < dimension; d++)
            position[d] = urand(minPosition[d], maxPosition[d]);
        Real minMomentum = -10;
        Real maxMomentum = 10;
        MomentumType momentum(urand(minMomentum, maxMomentum),
            urand(minMomentum, maxMomentum), urand(minMomentum, maxMomentum));
        FactorType factor = static_cast<FactorType>(urand(1e-5, 1e5));
        return Particle(position, momentum, factor, 0);
    }

    template<class ConstParticleRef1, class ConstParticleRef2>
    bool eqParticles_(ConstParticleRef1 a, ConstParticleRef2 b) const
    {
        return (a.getPosition() == b.getPosition()) &&
            (a.getMomentum() == b.getMomentum()) &&
            (a.getMass() == b.getMass()) &&
            (a.getCharge() == b.getCharge()) &&
            (a.getFactor() == b.getFactor());
    }
};


#endif
