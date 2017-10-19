#include "TestingUtility.h"

#include "pica/math/Constants.h"
#include "pica/particles/Particle.h"

#include <iostream>

using namespace pica;
using namespace std;


class ParticleTest: public BaseParticleFixture {

};

template <class ParticleType_>
class ParticleTest_ : public BaseFixture {
public:
    typedef ParticleType_ Particle;
    typedef typename Particle::PositionType PositionType;
    typedef typename Particle::MomentumType MomentumType;

    int getDimension() const
    {
        return VectorDimensionHelper<PositionType>::dimension;;
    }

    int getMomentumDimension() const
    {
        return VectorDimensionHelper<MomentumType>::dimension;;
    }

};


typedef ::testing::Types</*Particle1d,*/ Particle2d, Particle3d> types;
TYPED_TEST_CASE(ParticleTest_, types);


TYPED_TEST(ParticleTest_, DefaultConstructor)
{
    Particle particle;
    particle.setMass(Constants<double>::electronMass());
    ASSERT_EQ_VECTOR(PositionType(), particle.getPosition(), getDimension());
    ASSERT_EQ_VECTOR(MomentumType(), particle.getMomentum(), getMomentumDimension());
    ASSERT_EQ_VECTOR(MomentumType(), particle.getVelocity(), getMomentumDimension());
    ASSERT_EQ(1.0f, particle.getFactor());
    ASSERT_EQ(1.0, particle.getGamma());
}

TEST_F(ParticleTest, Constructor) {
    FP3 position(3.1, -32.1, 4.3e-5);
    FP3 momentum(-231.3e9, 0.0, 1.23e-5);
    int typeIndex = 1;
    float factor = 1.4e2f;
    Particle particle(position, momentum, typeIndex, factor);
    ASSERT_EQ_FP3(position, particle.getPosition());
    ASSERT_EQ_FP3(momentum, particle.getMomentum());
    ASSERT_EQ(typeIndex, particle.getType());
    ASSERT_EQ(factor, particle.getFactor());
    double expectedGamma = sqrt((FP)1 + momentum.norm2() / sqr(particle.mass() * constants::c));
    maxRelativeError = 1e-12;
    ASSERT_NEAR_FP(expectedGamma, particle.gamma());
}

TEST_F(ParticleTest, ConstructorDefaultFactor) {
    FP3 position(-12.34, 0.2, 423.12e-2);
    FP3 momentum(3254.23, -123.324, 1.23e5);
    int typeIndex = 0;
    Particle particle(position, momentum, typeIndex);
    ASSERT_EQ_FP3(position, particle.getPosition());
    ASSERT_EQ_FP3(momentum, particle.getMomentum());
    ASSERT_EQ(typeIndex, particle.getType());
    ASSERT_EQ(1.0f, particle.getFactor());
}

TEST_F(ParticleTest, CopyConstructor) {
    FP3 position(-134.12, 412.6342, 2346.562);
    FP3 momentum(-4531.23e5, 6534.123e3, 12.32);
    int typeIndex = 1;
    float factor = 213.51f;
    Particle particle(position, momentum, typeIndex, factor);
    Particle copyParticle(particle);
    ASSERT_TRUE(eqParticles(particle, copyParticle));
}

TEST_F(ParticleTest, Assignment) {
    FP3 position(432.453, -3452.15, -15.125);
    FP3 momentum(431.124, -54.12, 5643.176);
    int typeIndex = 0;
    float factor = 1.9945f;
    Particle particle(position, momentum, typeIndex, factor);
    Particle copyParticle;
    copyParticle = particle;
    ASSERT_TRUE(eqParticles(particle, copyParticle));
}

TEST_F(ParticleTest, Mass) {
    Particle particle = randomParticle();
    ASSERT_EQ(Particle::types[particle.getType()].mass, particle.mass());
}

TEST_F(ParticleTest, Charge) {
    Particle particle = randomParticle();
    ASSERT_EQ(Particle::types[particle.getType()].charge, particle.charge());
}

TEST_F(ParticleTest, GetType) {
    for (int i = 0; i < Particle::numTypes; i++) {
        FP3 position(-34.34, -45.12e3, 43.65e-2);
        FP3 momentum(231.123, -231.12, -432.12);
        Particle particle(position, momentum, i);
        ASSERT_EQ(i, particle.getType());
    }
}

TEST_F(ParticleTest, SetType) {
    for (int i = 0; i < Particle::numTypes; i++) {
        Particle particle = randomParticle();
        int newType = i;
        particle.setType(newType);
        ASSERT_EQ(newType, particle.getType());
        ASSERT_EQ(Particle::types[newType].mass, particle.mass());
        ASSERT_EQ(Particle::types[newType].charge, particle.charge());
    }
}

TEST_F(ParticleTest, GetSetPosition) {
    Particle particle = randomParticle();
    FP3 newPosition(54.126, -431.35, 35.65);
    particle.setPosition(newPosition);
    ASSERT_EQ_FP3(newPosition, particle.getPosition());
}

TEST_F(ParticleTest, GetSetMomentum) {
    Particle particle = randomParticle();
    FP3 newMomentum(54.12e+4, -543.63e-2, 643.165e5);
    particle.setMomentum(newMomentum);
    maxRelativeError = 1e-12;
    ASSERT_NEAR_FP3(newMomentum, particle.getMomentum());
    double expectedGamma = sqrt((FP)1 + newMomentum.norm2() / sqr(particle.mass() * constants::c));
    ASSERT_NEAR_FP(expectedGamma, particle.gamma());
}

TEST_F(ParticleTest, GetSetVelocity) {
    Particle particle = randomParticle();
    FP3 newVelocity(5243.1654, -56.23e5, -65.237e-4);
    particle.setVelocity(newVelocity);
    maxRelativeError = 1e-12;
    ASSERT_NEAR_FP3(newVelocity, particle.getVelocity());
    ///double expectedGamma = sqrt((FP)1 + newMomentum.norm2() / sqr(particle.mass() * Constants::c));
    ///ASSERT_NEAR_FP(expectedGamma, particle.gamma());
}

TEST_F(ParticleTest, Gamma) {
    Particle particle = randomParticle();
    double expectedGamma = sqrt((FP)1 + particle.getMomentum().norm2() /
                                        sqr(particle.mass() * constants::c));
    maxRelativeError = 1e-12;
    ASSERT_NEAR_FP(expectedGamma, particle.gamma());
}

TEST_F(ParticleTest, GetSetFactor) {
    Particle particle = randomParticle();
    float newFactor = 6523.54f;
    particle.setFactor(newFactor);
    ASSERT_EQ(newFactor, particle.getFactor());
}


