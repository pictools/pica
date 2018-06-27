#include "TestingUtility.h"

#include "pica/math/Constants.h"
#include "pica/particles/Particle.h"

using namespace pica;

template <class ParticleType>
class ParticleTest : public BaseParticleFixture<ParticleType> {
public:
    using typename BaseParticleFixture<ParticleType>::Particle;
    using typename BaseParticleFixture<ParticleType>::PositionType;
    using typename BaseParticleFixture<ParticleType>::MomentumType;
    using typename BaseParticleFixture<ParticleType>::GammaType;
    using typename BaseParticleFixture<ParticleType>::FactorType;
};

typedef ::testing::Types<Particle1d, Particle2d, Particle3d> types;
TYPED_TEST_CASE(ParticleTest, types);

TYPED_TEST(ParticleTest, DefaultConstructor)
{
    typedef typename ParticleTest<TypeParam>::Particle ParticleType;
    typedef typename ParticleTest<TypeParam>::PositionType PositionType;
    typedef typename ParticleTest<TypeParam>::MomentumType MomentumType;
    typedef typename ParticleTest<TypeParam>::GammaType GammaType;
    typedef typename ParticleTest<TypeParam>::FactorType FactorType;

    ParticleType particle;
    ASSERT_EQ_VECTOR(PositionType(), particle.getPosition(), this->dimension);
    ASSERT_EQ_VECTOR(MomentumType(), particle.getMomentum(), this->momentumDimension);
    ASSERT_EQ_VECTOR(MomentumType(), particle.getVelocity(), this->momentumDimension);
    ASSERT_EQ(static_cast<FactorType>(1.0), particle.getFactor());
    ASSERT_EQ(static_cast<GammaType>(1.0), particle.getGamma());
}

TYPED_TEST(ParticleTest, Constructor)
{
    typedef typename ParticleTest<TypeParam>::Particle ParticleType;
    typedef typename ParticleTest<TypeParam>::PositionType PositionType;
    typedef typename ParticleTest<TypeParam>::MomentumType MomentumType;
    typedef typename ParticleTest<TypeParam>::GammaType GammaType;
    typedef typename ParticleTest<TypeParam>::FactorType FactorType;

    PositionType position = this->getPosition(3.1, -32.1, 4.3e-5);
    MomentumType momentum(-231.3e9, 0.0, 1.23e-5);
    MassType mass = Constants<MassType>::electronMass();
    ChargeType charge = Constants<ChargeType>::electronCharge();
    FactorType factor = static_cast<FactorType>(1.4e2);

    ParticleType particle(position, momentum, factor, 0);
    ASSERT_EQ_VECTOR(position, particle.getPosition(), this->dimension);
    ASSERT_EQ_VECTOR(momentum, particle.getMomentum(), this->momentumDimension);
    ASSERT_EQ(mass, particle.getMass());
    ASSERT_EQ(charge, particle.getCharge());
    ASSERT_EQ(factor, particle.getFactor());
    GammaType expectedGamma = sqrt((FP)1 + momentum.norm2() / sqr(mass * Constants<GammaType >::c()));
    this->maxRelativeError = 1e-12;
    ASSERT_NEAR_FP(expectedGamma, particle.getGamma());
}

TYPED_TEST(ParticleTest, ConstructorDefaultFactor)
{
    typedef typename ParticleTest<TypeParam>::Particle ParticleType;
    typedef typename ParticleTest<TypeParam>::PositionType PositionType;
    typedef typename ParticleTest<TypeParam>::MomentumType MomentumType;
    typedef typename ParticleTest<TypeParam>::GammaType GammaType;
    typedef typename ParticleTest<TypeParam>::FactorType FactorType;

    PositionType position = this->getPosition(-12.34, 0.2, 423.12e-2);
    MomentumType momentum(3254.23, -123.324, 1.23e5);
    MassType mass = Constants<MassType>::electronMass();
    ChargeType charge = Constants<ChargeType>::electronCharge();
    ParticleType particle(position, momentum);
    ASSERT_EQ_VECTOR(position, particle.getPosition(), this->dimension);
    ASSERT_EQ_VECTOR(momentum, particle.getMomentum(), this->momentumDimension);
    ASSERT_EQ(mass, particle.getMass());
    ASSERT_EQ(charge, particle.getCharge());
    ASSERT_EQ(static_cast<FactorType>(1.0), particle.getFactor());
}

TYPED_TEST(ParticleTest, CopyConstructor)
{
    typedef typename ParticleTest<TypeParam>::Particle ParticleType;
    typedef typename ParticleTest<TypeParam>::PositionType PositionType;
    typedef typename ParticleTest<TypeParam>::MomentumType MomentumType;
    typedef typename ParticleTest<TypeParam>::GammaType GammaType;
    typedef typename ParticleTest<TypeParam>::FactorType FactorType;

    PositionType position = this->getPosition(-134.12, 412.6342, 2346.562);
    MomentumType momentum(-4531.23e5, 6534.123e3, 12.32);
    FactorType factor = 213.51;
    ParticleType particle(position, momentum, factor, 0);
    ParticleType copyParticle(particle);
    ASSERT_TRUE(this->eqParticles_(particle, copyParticle));
}

TYPED_TEST(ParticleTest, Assignment)
{
    typedef typename ParticleTest<TypeParam>::Particle ParticleType;
    typedef typename ParticleTest<TypeParam>::PositionType PositionType;
    typedef typename ParticleTest<TypeParam>::MomentumType MomentumType;
    typedef typename ParticleTest<TypeParam>::GammaType GammaType;
    typedef typename ParticleTest<TypeParam>::FactorType FactorType;

    PositionType position = this->getPosition(432.453, -3452.15, -15.125);
    MomentumType momentum(431.124, -54.12, 5643.176);
    FactorType factor = 1.9945;
    ParticleType particle(position, momentum, factor, 0);
    ParticleType copyParticle;
    copyParticle = particle;
    ASSERT_TRUE(this->eqParticles_(particle, copyParticle));
}

TYPED_TEST(ParticleTest, GetSetPosition)
{
    typedef typename ParticleTest<TypeParam>::Particle ParticleType;
    typedef typename ParticleTest<TypeParam>::PositionType PositionType;

    ParticleType particle = this->randomParticle();
    PositionType newPosition = this->getPosition(54.126, -431.35, 35.65);
    particle.setPosition(newPosition);
    ASSERT_EQ_VECTOR(newPosition, particle.getPosition(), this->dimension);
}

TYPED_TEST(ParticleTest, GetSetMomentum)
{
    typedef typename ParticleTest<TypeParam>::Particle ParticleType;
    typedef typename ParticleTest<TypeParam>::MomentumType MomentumType;
    typedef typename ParticleTest<TypeParam>::GammaType GammaType;

    ParticleType particle = this->randomParticle();
    MomentumType newMomentum(54.12e+4, -543.63e-2, 643.165e5);
    particle.setMomentum(newMomentum);
    this->maxRelativeError = 1e-12;
    ASSERT_NEAR_VECTOR(newMomentum, particle.getMomentum());
    GammaType expectedGamma = sqrt((FP)1 + newMomentum.norm2() / sqr(particle.getMass() * Constants<GammaType >::c()));
    ASSERT_NEAR_FP(expectedGamma, particle.getGamma());
}

TYPED_TEST(ParticleTest, GetSetVelocity)
{
    typedef typename ParticleTest<TypeParam>::Particle ParticleType;
    typedef typename ParticleTest<TypeParam>::MomentumType MomentumType;

    ParticleType particle = this->randomParticle();
    MomentumType newVelocity(5243.1654, -56.23e5, -65.237e-4);
    particle.setVelocity(newVelocity);
    this->maxRelativeError = 1e-12;
    MomentumType v = particle.getVelocity();
    ASSERT_NEAR_VECTOR(newVelocity, particle.getVelocity());
}

TYPED_TEST(ParticleTest, GetGamma)
{
    typedef typename ParticleTest<TypeParam>::Particle ParticleType;
    typedef typename ParticleTest<TypeParam>::GammaType GammaType;

    ParticleType particle = this->randomParticle();
    GammaType expectedGamma = sqrt(static_cast<GammaType>(1.0) + particle.getMomentum().norm2() /
        sqr(particle.getMass() * Constants<GammaType>::c()));
    this->maxRelativeError = 1e-12;
    ASSERT_NEAR_FP(expectedGamma, particle.getGamma());
}

TYPED_TEST(ParticleTest, GetMass)
{
    typedef typename ParticleTest<TypeParam>::Particle ParticleType;

    ParticleType particle = this->randomParticle();
    ASSERT_EQ(pica::Constants<MassType>::electronMass(), particle.getMass());
}

TYPED_TEST(ParticleTest, GetCharge)
{
    typedef typename ParticleTest<TypeParam>::Particle ParticleType;

    ParticleType particle = this->randomParticle();
    ASSERT_EQ(pica::Constants<ChargeType>::electronCharge(), particle.getCharge());
}

TYPED_TEST(ParticleTest, GetSetFactor) 
{
    typedef typename ParticleTest<TypeParam>::Particle ParticleType;
    typedef typename ParticleTest<TypeParam>::FactorType FactorType;

    ParticleType particle = this->randomParticle();
    FactorType newFactor = 6523.54;
    particle.setFactor(newFactor);
    ASSERT_EQ(newFactor, particle.getFactor());
}
