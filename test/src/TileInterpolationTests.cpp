#include "TestingUtility.h"

#include "pica/Ensemble.h"
#include "pica/grid/Grid.h"
#include "pica/TileInterpolator.h"

using namespace pica;


class TileInterpolationTest : public BaseParticleFixture {
public:

    Grid* grid;
    Ensemble* ensemble;


    virtual void SetUp() {
        BaseParticleFixture::SetUp();
        maxAbsoluteError = (FP)1e-4;
        maxRelativeError = (FP)0.001;
        initGrid();
        initEnsemble();       
    }

    void initGrid() {
        grid = new Grid(parameters.gridSize, parameters.timeStep,
            parameters.localMin, parameters.step, parameters.gridSize, Int3(0, 0, 0));
        for (int i = 0; i < grid->numCells.x; i++)
        for (int j = 0; j < grid->numCells.y; j++)
        for (int k = 0; k < grid->numCells.z; k++) {
            grid->Ex(i, j, k) = grid->ExPosition(i, j, k).x;
            grid->Ey(i, j, k) = grid->EyPosition(i, j, k).y;
            grid->Ez(i, j, k) = grid->EzPosition(i, j, k).z;
            grid->Bx(i, j, k) = grid->BxPosition(i, j, k).x;
            grid->By(i, j, k) = grid->ByPosition(i, j, k).y;
            grid->Bz(i, j, k) = grid->BzPosition(i, j, k).z;
        }
    }

    void initEnsemble() {
        ensemble = new Ensemble(parameters.gridSize, parameters.localMin, parameters.step, 3);
        for (int i = 0; i < 100; i++)
            ensemble->add(randomParticle());
    }

    virtual void TearDown() {
        delete grid;
        delete ensemble;
    }

    template<typename TileInterpolator>
    void doTest() {
        maxAbsoluteError = (FP)1e-10;
        maxRelativeError = (FP)1e-10;
        for (int i = 1; i < ensemble->numSystems.x - 1; i++)
        for (int j = 1; j < ensemble->numSystems.y - 1; j++)
        for (int k = 1; k < ensemble->numSystems.z - 1; k++) {
            TileInterpolator interpolator(i, j, k, grid);
            ParticleSystem* system = ensemble->system(i, j, k);
            for (system->begin(); !system->end(); system->next()) {
                FP3 coords = system->currentParticle->getPosition();
                FP3 e, b, tileE, tileB;
                grid->getFields(coords.x, coords.y, coords.z, e, b);
                interpolator.getFields(coords, tileE, tileB);
                ASSERT_NEAR_FP3(e, tileE);
                ASSERT_NEAR_FP3(b, tileB);
            }
        }
    }

};


TEST_F(TileInterpolationTest, CIC) {
    grid->setInterpolationType(Grid::Interpolation_CIC);
    doTest<TileInterpolatorCIC>();
}


TEST_F(TileInterpolationTest, TSC) {
    grid->setInterpolationType(Grid::Interpolation_TSC);
    doTest<TileInterpolatorTSC>();
}
