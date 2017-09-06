#include "TestingUtility.h"

#include "pica/math/FP.h"

using namespace pica;


TEST(FPTest, Sqr) {
    FP args[] = {3.1, -2.41, 1.1, 0.0, -7.21};
    for (int i = 0; i < sizeof(args) / sizeof(*args); i++)
        ASSERT_EQ(sqr(args[i]), args[i] * args[i]);
}
