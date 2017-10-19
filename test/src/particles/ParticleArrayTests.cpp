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
            if (!eqParticles(a[i], b[i]))
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

TYPED_TEST(ParticleArrayTest, PushBack)
{
    ParticleArray particles;
    for (int i = 0; i < 17; i++)
        particles.pushBack(randomParticle());
    ParticleArray particlesCopy(particles);
    ASSERT_TRUE(eqParticleArrays(particles, particlesCopy)); 
}
