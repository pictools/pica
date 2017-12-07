#include "TestingUtility.h"

#include "pica/math/Dimension.h"

using namespace pica;


TEST(DimensionTest, Values) {
    ASSERT_EQ(One, 1);
    ASSERT_EQ(Two, 2);
    ASSERT_EQ(Three, 3);
}
