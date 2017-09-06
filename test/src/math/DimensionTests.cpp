#include "TestingUtility.h"

#include "pica/math/Dimension.h"

using namespace pica;


TEST(DimensionTest, Values) {
    ASSERT_EQ(Dimension::One, 1);
    ASSERT_EQ(Dimension::Two, 2);
    ASSERT_EQ(Dimension::Three, 3);
}
