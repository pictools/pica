#include "TestingUtility.h"

#include "pica/fieldInterpolation/FieldInterpolator.h"
#include "pica/grid/YeeGrid.h"

using namespace pica;


TEST(FieldInterpolatorTest, Constructor)
{
    double origin = -1.5;
    double step = 0.5;
    int size = 10;
    YeeGrid<One, double> grid(origin, step, size);
    FieldInterpolatorCIC<YeeGrid<One, double> > interpolator(grid);
    Vector3<double> e, b;
    interpolator.get(1.0, e, b);

    YeeGrid<Two, double> grid2(Vector2<double>(1.0, 1.0), Vector2<double>(2.0, 2.0), Vector2<int>(3, 5));
    FieldInterpolatorCIC<YeeGrid<Two, double> > interpolator2(grid2);

    YeeGrid<Three, double> grid3(Vector3<double>(1.0, 1.0, 1.0), Vector3<double>(2.0, 2.0, 2.0), Vector3<int>(3, 5, 7));
    FieldInterpolatorCIC<YeeGrid<Three, double> > interpolator3(grid3);
}