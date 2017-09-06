#include "TestingUtility.h"

#include "pica/math/Constants.h"
#include "pica/grid/Grid.h"
#include "pica/particles/Ensemble.h"
#include "pica/utility/Utility.h"

#include <iostream>
using std::cout;

using namespace pica;


class ParticleMoveTest : public BaseParticleFixture
{

protected:

    virtual void SetUp();

};


void ParticleMoveTest::SetUp()
{
    BaseParticleFixture::SetUp();
    maxAbsoluteError = (FP)1e-12;
    maxRelativeError = (FP)0.01;

    // in this group of tests we need gigantic physical space
    // as consequently steps are also gigantic - it's OK for these tests
    Parameters parameters;
    FP spaceSize = (FP)1e10;
    parameters.localMin = parameters.globalMin = FP3(-1, -1, -1) * spaceSize;
    parameters.localMax = parameters.globalMax = FP3(1, 1, 1) * spaceSize;
    parameters.localGridSize = parameters.globalGridSize = parameters.gridSize = Int3(32, 32, 32);
    parameters.step = (parameters.globalMax - parameters.globalMin) /
        FP3(parameters.globalGridSize);
    parameters.localMinIndex = Int3(0, 0, 0);
    parameters.localMaxIndex = parameters.localGridSize - Int3(1, 1, 1);
    parameters.numIterations = 1;
    parameters.timeStep = 0;
}
