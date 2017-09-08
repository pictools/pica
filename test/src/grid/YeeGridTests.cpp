#include "TestingUtility.h"

#include "pica/grid/YeeGrid.h"

using namespace pica;


TEST(YeeGridTest, Constructor)
{
    YeeGrid<One, double> grid1d(3);
    YeeGrid<Two, float> grid2d(YeeGrid<Two, float>::IndexType(5, 6));
    YeeGrid<Three> grid3d(YeeGrid<Three>::IndexType(4, 3, 9));
}

