#include "TestingUtility.h"

#include "pica/particles/Particle.h"
#include "pica/particles/ParticleArray.h"


using namespace pica;


template <class ParticleArrayType>
class ParticleArrayTest : public BaseParticleFixture_<typename ParticleArrayType::Particle> {
public:
    typedef ParticleArrayType ParticleArray;

    bool eqParticleArrays(ParticleArray& a, ParticleArray& b) const
    {
        if (a.size() != b.size())
            return false;
        for (int i = 0; i < a.size(); i++)
            if (!eqParticles_(a[i], b[i]))
                return false;
        return true;
    }

};

typedef ::testing::Types<
    ParticleArray<Particle1d, ParticleRepresentation_AoS>::Type,
    ParticleArray<Particle2d, ParticleRepresentation_AoS>::Type,
    ParticleArray<Particle3d, ParticleRepresentation_AoS>::Type,
    ParticleArray<Particle1d, ParticleRepresentation_SoA>::Type,
    ParticleArray<Particle2d, ParticleRepresentation_SoA>::Type,
    ParticleArray<Particle3d, ParticleRepresentation_SoA>::Type
> types;
TYPED_TEST_CASE(ParticleArrayTest, types);

TYPED_TEST(ParticleArrayTest, DefaultConstructor)
{
    ParticleArray particles;
    ASSERT_EQ(0, particles.size());
}

TYPED_TEST(ParticleArrayTest, Assignment)
{
    ParticleArray particles, particlesCopy, particlesAnotherCopy;
    for (int i = 0; i < 9; i++)
        particles.pushBack(randomParticle());
    for (int i = 0; i < 11; i++)
        particlesCopy.pushBack(randomParticle());
    particlesCopy = particles;
    particlesAnotherCopy = particles;
    ASSERT_TRUE(eqParticleArrays(particles, particlesCopy));
    ASSERT_TRUE(eqParticleArrays(particles, particlesAnotherCopy));
}

TYPED_TEST(ParticleArrayTest, Size)
{
    ParticleArray particles;
    int numParticles = 12;
    for (int i = 0; i < numParticles; i++)
        particles.pushBack(randomParticle());
    ASSERT_EQ(numParticles, particles.size());
}

TYPED_TEST(ParticleArrayTest, IndexAccess)
{
    ParticleArray particles;
    const int numParticles = 15;
    Particle particleArray[numParticles];
    for (int i = 0; i < numParticles; i++) {
        Particle particle = randomParticle();
        particleArray[i] = particle;
        particles.pushBack(particle);
    }
    const ParticleArray constParticles = particles;
    ASSERT_EQ(numParticles, particles.size());
    ASSERT_EQ(numParticles, constParticles.size());
    for (int i = 0; i < numParticles; i++) {
        EXPECT_TRUE(eqParticles_(particleArray[i], particles[i]));
        EXPECT_TRUE(eqParticles_(particleArray[i], constParticles[i]));
    }
}

TYPED_TEST(ParticleArrayTest, Back)
{
    ParticleArray particles;
    for (int i = 0; i < 13; i++) {
        Particle particle = randomParticle();
        particles.pushBack(particle);
        EXPECT_TRUE(eqParticles_(particle, particles.back()));
        const ParticleArray particlesCopy = particles;
        EXPECT_TRUE(eqParticles_(particle, particlesCopy.back()));
    }
}

TYPED_TEST(ParticleArrayTest, PushBack)
{
    ParticleArray particles;
    for (int i = 0; i < 17; i++)
        particles.pushBack(randomParticle());
    ParticleArray particlesCopy(particles);
    ASSERT_TRUE(eqParticleArrays(particles, particlesCopy));
}

TYPED_TEST(ParticleArrayTest, PopBack)
{
    ParticleArray particles;
    for (int i = 0; i < 33; i++)
        particles.pushBack(randomParticle());
    for (int i = particles.size(); i > 0; i--) {
        ParticleArray particlesCopy = particles;
        particles.popBack();
        ASSERT_EQ(particlesCopy.size(), particles.size() + 1);
        for (int j = 0; j < particles.size(); j++)
            EXPECT_TRUE(eqParticles_(particles[j], particlesCopy[j]));
    }
    ParticleArray particlesCopy(particles);
    ASSERT_TRUE(eqParticleArrays(particles, particlesCopy));
}
