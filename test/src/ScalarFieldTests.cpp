#include "TestingUtility.h"

#include "pica/ScalarField.h"

using namespace pica;


class ScalarFieldTest : public BaseFixture {
protected:
    virtual void SetUp() {
        BaseFixture::SetUp();
        maxAbsoluteError = (FP)1e-4;
        maxRelativeError = (FP)0.001;
    }
};

TEST_F(ScalarFieldTest, Constructor) {
    ScalarField f(Int3(3, 4, 5));
    f(0, 0, 0) = 1;
    f(2, 3, 4) = 2;
}

TEST_F(ScalarFieldTest, CopyConstructor) {
    ScalarField f(Int3(3, 4, 5));
    ScalarField g(f);
    g(0, 0, 0) = 1;
    g(2, 3, 4) = 2;
}

TEST_F(ScalarFieldTest, Assignment) {
    ScalarField f(Int3(3, 4, 5)), g(Int3(1, 3, 2));
    g = f;
    g(0, 0, 0) = 1;
    g(2, 3, 4) = 2;
}

ScalarField createScalarField(const Int3& size) {
    ScalarField f(size);
    for (int i = 0; i < size.x; i++)
    for (int j = 0; j < size.y; j++)
    for (int k = 0; k < size.z; k++)
        f(i, j, k) = k + (j + i * size.y) * size.z;
    return f;
}

TEST_F(ScalarFieldTest, IndexAccess) {
    Int3 size(5, 3, 8);
    ScalarField f1(createScalarField(size));
    const ScalarField f2(createScalarField(size));
    for (int i = 0; i < size.x; i++)
    for (int j = 0; j < size.y; j++)
    for (int k = 0; k < size.z; k++) {
        ASSERT_EQ(f1(i, j, k), k + (j + i * size.y) * size.z);
        ASSERT_EQ(f1(i, j, k), f2(i, j, k));
        ASSERT_EQ(f1(i, j, k), f1(Int3(i, j, k)));
        ASSERT_EQ(f1(i, j, k), f2(Int3(i, j, k)));
    }
}

TEST_F(ScalarFieldTest, depositCIC) {
    Int3 size(5, 8, 6);
    ScalarField f(size);
    FP old = 0;
    for (int testIdx = 0; testIdx < 100; testIdx++) {
        Int3 base = urandInt3(Int3(0, 0, 0), size - Int3(2, 2, 2));
        FP3 coeffs = urandFP3(FP3(0, 0, 0), FP3(1, 1, 1));
        FP value = urand(1, 10);
        f.depositCIC(value, base, coeffs);
        FP actual = 0;
        for(int i = 0; i < size.x; i++)
        for(int j = 0; j < size.y; j++)
        for(int k = 0; k < size.z; k++)
            actual += f(i, j, k);
        FP diff = actual - old;
        ASSERT_NEAR_FP(value, diff);
        old = actual;
    }
}

TEST_F(ScalarFieldTest, depositTSC) {
    Int3 size(5, 8, 6);
    ScalarField f(size);
    FP old = 0;
    for (int testIdx = 0; testIdx < 100; testIdx++) {
        Int3 base = urandInt3(Int3(1, 1, 1), size - Int3(2, 2, 2));
        FP3 coeffs = urandFP3(FP3(0, 0, 0), FP3(1, 1, 1));
        FP value = urand(1, 10);
        f.depositTSC(value, base, coeffs);
        FP actual = 0;
        for(int i = 0; i < size.x; i++)
        for(int j = 0; j < size.y; j++)
        for(int k = 0; k < size.z; k++)
            actual += f(i, j, k);
        FP diff = actual - old;
        ASSERT_NEAR_FP(value, diff);
        old = actual;
    }
}

TEST_F(ScalarFieldTest, depositPCS) {
    Int3 size(7, 10, 9);
    ScalarField f(size);
    FP old = 0;
    for (int testIdx = 0; testIdx < 100; testIdx++) {
        Int3 base = urandInt3(Int3(1, 1, 1), size - Int3(3, 3, 3));
        FP3 coeffs = urandFP3(FP3(0, 0, 0), FP3(1, 1, 1));
        FP value = urand(1, 10);
        f.depositPCS(value, base, coeffs);
        FP actual = 0;
        for(int i = 0; i < size.x; i++)
        for(int j = 0; j < size.y; j++)
        for(int k = 0; k < size.z; k++)
            actual += f(i, j, k);
        FP diff = actual - old;
        ASSERT_NEAR_FP(value, diff);
        old = actual;
    }
}

TEST_F(ScalarFieldTest, Zeroize) {
    Int3 size(5, 3, 8);
    ScalarField f(createScalarField(size));
    f(0, 0, 0) = 1;
    f.zeroize();
    for (int i = 0; i < size.x; i++)
    for (int j = 0; j < size.y; j++)
    for (int k = 0; k < size.z; k++)
        ASSERT_EQ(f(i, j, k), 0);
}
