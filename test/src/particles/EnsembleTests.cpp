#include "TestingUtility.h"

#include "pica/particles/Ensemble.h"
#include "pica/particles/ParticleArray.h"
#include "pica/utility/Utility.h"

#include <string>
#include <vector>

using namespace pica;
using std::string;
using std::vector;


template <class EnsembleType>
class EnsembleTest : public BaseParticleFixture_<typename EnsembleType::Particle> {
public:
    typedef EnsembleType Ensemble;
};

typedef ::testing::Types<
    EnsembleUnordered<ParticleArray<Particle1d, ParticleRepresentation_AoS>::Type>,
    EnsembleUnordered<ParticleArray<Particle2d, ParticleRepresentation_AoS>::Type>,
    EnsembleUnordered<ParticleArray<Particle3d, ParticleRepresentation_AoS>::Type>,
    EnsembleUnordered<ParticleArray<Particle1d, ParticleRepresentation_SoA>::Type>,
    EnsembleUnordered<ParticleArray<Particle2d, ParticleRepresentation_SoA>::Type>,
    EnsembleUnordered<ParticleArray<Particle3d, ParticleRepresentation_SoA>::Type>,
    EnsembleOrdered<ParticleArray<Particle1d, ParticleRepresentation_AoS>::Type>,
    EnsembleOrdered<ParticleArray<Particle2d, ParticleRepresentation_AoS>::Type>,
    EnsembleOrdered<ParticleArray<Particle3d, ParticleRepresentation_AoS>::Type>,
    EnsembleOrdered<ParticleArray<Particle1d, ParticleRepresentation_SoA>::Type>,
    EnsembleOrdered<ParticleArray<Particle2d, ParticleRepresentation_SoA>::Type>,
    EnsembleOrdered<ParticleArray<Particle3d, ParticleRepresentation_SoA>::Type>
> types;
TYPED_TEST_CASE(EnsembleTest, types);

TYPED_TEST(EnsembleTest, DefaultConstructor)
{
    Ensemble ensemble(getPosition(0, 0, 0), getPosition(1, 1, 1));
    ASSERT_EQ(0, ensemble.size());
}

TYPED_TEST(EnsembleTest, Add)
{
    PositionType minPosition = getPosition(0, 0, 0);
    PositionType maxPosition = getPosition(1, 1, 1);
    Ensemble ensemble(minPosition, maxPosition);
    int numParticles = 37;
    for (int i = 0; i < numParticles; i++)
        ensemble.add(randomParticle(minPosition, maxPosition));
    EXPECT_EQ(numParticles, ensemble.size());
}
