#include "TestingUtility.h"

#include "pica/fieldSolver/YeeSolver.h"
#include "pica/grid/YeeGrid.h"
#include "pica/math/Vectors.h"

using namespace pica;


TEST(YeeSolverTest, Constructor)
{
    YeeSolver solver;
    solver;
}

TEST(YeeSolverTest, UpdateE1d)
{
    YeeSolver solver;
    YeeGrid<One> grid(Vector1<double>(0.0), Vector1<double>(0.1), Vector1<int>(10));
    solver.updateE(grid, 0.1);
}

TEST(YeeSolverTest, UpdateB1d)
{
    YeeSolver solver;
    YeeGrid<One> grid(Vector1<double>(0.0), Vector1<double>(0.1), Vector1<int>(10));
    solver.updateB(grid, 0.1);
}

TEST(YeeSolverTest, UpdateE2d)
{
    YeeSolver solver;
    YeeGrid<Two, float> grid(Vector2<float>(0.0f, 0.0f), Vector2<float>(0.1f, 0.1f), Vector2<int>(5, 7));
    solver.updateE(grid, 0.1f);
}

TEST(YeeSolverTest, UpdateB2d)
{
    YeeSolver solver;
    YeeGrid<Two, float> grid(Vector2<float>(0.0f, 0.0f), Vector2<float>(0.1f, 0.1f), Vector2<int>(5, 7));
    solver.updateB(grid, 0.1f);
}

TEST(YeeSolverTest, UpdateE3d)
{
    YeeSolver solver;
    YeeGrid<Three, double> grid(Vector3<double>(0.0, 0.0, 0.0), Vector3<double>(0.1, 0.1, 0.1), Vector3<int>(5, 7, 10));
    solver.updateE(grid, 0.2);
}

TEST(YeeSolverTest, UpdateB3d)
{
    YeeSolver solver;
    YeeGrid<Three, double> grid(Vector3<double>(0.0, 0.0, 0.0), Vector3<double>(0.1, 0.1, 0.1), Vector3<int>(5, 7, 10));
    solver.updateB(grid, 0.2);
}
