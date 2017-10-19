#include "TestingUtility.h"

#include "pica/math/Constants.h"
#include "pica/particles/Particle.h"
#include "pica/particles/ParticleTraits.h"

#include <iostream>

using namespace pica;
using namespace std;


template <class ParticleType>
class ParticleTest : public BaseFixture {
public:
    typedef ParticleType Particle;
    typedef typename ParticleTraits<Particle>::PositionType PositionType;
    typedef typename ParticleTraits<Particle>::MomentumType MomentumType;
    typedef typename ParticleTraits<Particle>::GammaType GammaType;
    typedef typename ParticleTraits<Particle>::MassType MassType;
    typedef typename ParticleTraits<Particle>::ChargeType ChargeType;
    typedef typename ParticleTraits<Particle>::FactorType FactorType;
    typedef typename ScalarType<MomentumType>::Type Real;

    int getDimension() const
    {
        return VectorDimensionHelper<PositionType>::dimension;;
    }

    int getMomentumDimension() const
    {
        return VectorDimensionHelper<MomentumType>::dimension;;
    }

    // Helper function to unify initialization of positions for 1d, 2d and 3d
    // In 1d y, z are ignored, in 2d z is ignored
    PositionType getPosition(Real x, Real y, Real z) const
    {
        Real positionArray[] = { x, y, z };
        PositionType position;
        for (int d = 0; d < getDimension(); d++)
            position[d] = positionArray[d];
        return position;
    }

    Particle randomParticle() const
    {
        Real minPosition = -10;
        Real maxPosition = 10;
        PositionType position = getPosition(urand(minPosition, maxPosition),
            urand(minPosition, maxPosition), urand(minPosition, maxPosition));
        Real minMomentum = -10;
        Real maxMomentum = 10;
        MomentumType momentum(urand(minMomentum, maxMomentum),
            urand(minMomentum, maxMomentum), urand(minMomentum, maxMomentum));
        FactorType factor = urand(1e-5, 1e5);
        return Particle(position, momentum, Constants<MassType>::electronMass(),
            Constants<ChargeType>::electronCharge(), factor);
    }

};

typedef ::testing::Types</*Particle1d,*/ Particle2d, Particle3d> types;
TYPED_TEST_CASE(ParticleTest, types);

TYPED_TEST(ParticleTest, DefaultConstructor)
{
    Particle particle;
    particle.setMass(Constants<double>::electronMass());
    ASSERT_EQ_VECTOR(PositionType(), particle.getPosition(), getDimension());
    ASSERT_EQ_VECTOR(MomentumType(), particle.getMomentum(), getMomentumDimension());
    ASSERT_EQ_VECTOR(MomentumType(), particle.getVelocity(), getMomentumDimension());
    ASSERT_EQ(static_cast<FactorType>(1.0), particle.getFactor());
    ASSERT_EQ(static_cast<GammaType>(1.0), particle.getGamma());
}

TYPED_TEST(ParticleTest, Constructor)
{
    PositionType position = getPosition(3.1, -32.1, 4.3e-5);
    MomentumType momentum(-231.3e9, 0.0, 1.23e-5);
    MassType mass = Constants<MassType>::electronMass();
    ChargeType charge = Constants<ChargeType>::electronCharge();
    FactorType factor = 1.4e2f;

    Particle particle(position, momentum, mass, charge, factor);
    ASSERT_EQ_VECTOR(position, particle.getPosition(), getDimension());
    ASSERT_EQ_VECTOR(momentum, particle.getMomentum(), getMomentumDimension());
    ASSERT_EQ(mass, particle.getMass());
    ASSERT_EQ(charge, particle.getCharge());
    ASSERT_EQ(factor, particle.getFactor());
    Real expectedGamma = sqrt((FP)1 + momentum.norm2() / sqr(mass * Constants<double>::c()));
    maxRelativeError = 1e-12;
    ASSERT_NEAR_FP(expectedGamma, particle.getGamma());
}

TYPED_TEST(ParticleTest, ConstructorDefaultFactor)
{
    PositionType position = getPosition(-12.34, 0.2, 423.12e-2);
    MomentumType momentum(3254.23, -123.324, 1.23e5);
    MassType mass = Constants<MassType>::electronMass();
    ChargeType charge = Constants<ChargeType>::electronCharge();
    Particle particle(position, momentum, mass, charge);
    ASSERT_EQ_VECTOR(position, particle.getPosition(), getDimension());
    ASSERT_EQ_VECTOR(momentum, particle.getMomentum(), getMomentumDimension());
    ASSERT_EQ(mass, particle.getMass());
    ASSERT_EQ(charge, particle.getCharge());
    ASSERT_EQ(static_cast<FactorType>(1.0), particle.getFactor());
}

TYPED_TEST(ParticleTest, CopyConstructor)
{
    PositionType position = getPosition(-134.12, 412.6342, 2346.562);
    MomentumType momentum(-4531.23e5, 6534.123e3, 12.32);
    MassType mass = Constants<MassType>::electronMass();
    ChargeType charge = Constants<ChargeType>::electronCharge();
    FactorType factor = 213.51f;
    Particle particle(position, momentum, mass, charge, factor);
    Particle copyParticle(particle);
    ASSERT_TRUE(eqParticles(particle, copyParticle));
}

TYPED_TEST(ParticleTest, Assignment)
{
    PositionType position = getPosition(432.453, -3452.15, -15.125);
    MomentumType momentum(431.124, -54.12, 5643.176);
    MassType mass = Constants<MassType>::electronMass();
    ChargeType charge = Constants<ChargeType>::electronCharge();
    FactorType factor = 1.9945f;
    Particle particle(position, momentum, mass, charge, factor);
    Particle copyParticle;
    copyParticle = particle;
    ASSERT_TRUE(eqParticles(particle, copyParticle));
}
 
TYPED_TEST(ParticleTest, GetSetPosition)
{
    Particle particle = randomParticle();
    PositionType newPosition = getPosition(54.126, -431.35, 35.65);
    particle.setPosition(newPosition);
    ASSERT_EQ_VECTOR(newPosition, particle.getPosition(), getDimension());
}

TYPED_TEST(ParticleTest, GetSetMomentum)
{
    Particle particle = randomParticle();
    MomentumType newMomentum(54.12e+4, -543.63e-2, 643.165e5);
    particle.setMomentum(newMomentum);
    maxRelativeError = 1e-12;
    ASSERT_NEAR_VECTOR(newMomentum, particle.getMomentum());
    double expectedGamma = sqrt((FP)1 + newMomentum.norm2() / sqr(particle.getMass() * constants::c));
    ASSERT_NEAR_FP(expectedGamma, particle.getGamma());
}

TYPED_TEST(ParticleTest, GetSetVelocity)
{
    Particle particle = randomParticle();
    MomentumType newVelocity(5243.1654, -56.23e5, -65.237e-4);
    particle.setVelocity(newVelocity);
    maxRelativeError = 1e-12;
    MomentumType v = particle.getVelocity();
    ASSERT_NEAR_VECTOR(newVelocity, particle.getVelocity());
}

TYPED_TEST(ParticleTest, GetGamma)
{
    Particle particle = randomParticle();
    GammaType expectedGamma = sqrt(static_cast<GammaType>(1.0) + particle.getMomentum().norm2() /
        sqr(particle.getMass() * Constants<GammaType>::c()));
    maxRelativeError = 1e-12;
    ASSERT_NEAR_FP(expectedGamma, particle.getGamma());
}

TYPED_TEST(ParticleTest, GetSetMass)
{
    Particle particle = randomParticle();
    MassType newMass = 1.3e-10;
    particle.setMass(newMass);
    ASSERT_EQ(newMass, particle.getMass());
}

TYPED_TEST(ParticleTest, GetSetCharge)
{
    Particle particle = randomParticle();
    ChargeType newCharge = -5.7e-13;
    particle.setCharge(newCharge);
    ASSERT_EQ(newCharge, particle.getCharge());
}

TYPED_TEST(ParticleTest, GetSetFactor) 
{
    Particle particle = randomParticle();
    FactorType newFactor = 6523.54f;
    particle.setFactor(newFactor);
    ASSERT_EQ(newFactor, particle.getFactor());
}
