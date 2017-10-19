#include "TestingUtility.h"

#include "pica/utility/Array.h"

#include <typeinfo>

using namespace pica;


TEST(Array1dTest, ValueType)
{
    ASSERT_TRUE(typeid(Array1d<float>::ValueType) == typeid(float));
}

TEST(Array1dTest, DefaultConstructor)
{
    Array1d<double> a;
    ASSERT_EQ(0, a.getSize().x);
}

TEST(Array1dTest, Constructor)
{
    Array1d<int> a(5);
    ASSERT_EQ(5, a.getSize().x);
    for (int i = 0; i < a.getSize().x; i++)
        ASSERT_EQ(0, a(i));
    Array1d<float> b(9, -3.0f);
    ASSERT_EQ(9, b.getSize().x);
    for (int i = 0; i < b.getSize().x; i++)
        ASSERT_EQ(-3.0f, b(i));
}

TEST(Array1dTest, GetSize)
{
    const Array1d<double> a(6);
    ASSERT_EQ(6, a.getSize().x); 
}

TEST(Array1dTest, Access)
{
    Array1d<int> a(4, -1);
    for (int i = 0; i < 4; i++)
        a(i) = i;
    for (int i = 0; i < 4; i++)
        ASSERT_EQ(i, a(i));
}

TEST(Array1dTest, AccessConst)
{
    const Array1d<int> a(5, 2);
    for (int i = 0; i < 5; i++)
        ASSERT_EQ(2, a(i));
}


TEST(Array2dTest, ValueType)
{
    ASSERT_TRUE(typeid(Array2d<int>::ValueType) == typeid(int));
}

TEST(Array2dTest, DefaultConstructor)
{
    Array2d<float> a;
    ASSERT_TRUE(Array2d<double>::IndexType(0, 0) == a.getSize());
}

TEST(Array2dTest, Constructor)
{
    Array2d<double> a(4, 3);
    ASSERT_TRUE(Array2d<double>::IndexType(4, 3) == a.getSize());
    for (int i = 0; i < a.getSize().x; i++)
    for (int j = 0; j < a.getSize().y; j++)
        ASSERT_EQ(0, a(i, j));
    Array2d<int> b(Array2d<int>::IndexType(5, 2), 3);
    ASSERT_TRUE(Array2d<double>::IndexType(5, 2) == b.getSize());
    for (int i = 0; i < b.getSize().x; i++)
    for (int j = 0; j < b.getSize().y; j++)
        ASSERT_EQ(3, b(i, j));
}

TEST(Array2dTest, GetSize)
{
    const Array2d<float> a(7, 4);
    ASSERT_TRUE(Array2d<float>::IndexType(7, 4) == a.getSize());
}

TEST(Array2dTest, Access)
{
    Array2d<int> a(3, 5, -1);
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 5; j++)
        a(i, j) = i + j;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 5; j++)
        ASSERT_EQ(i + j, a(i, j));
}

TEST(Array2dTest, AccessConst)
{
    const Array2d<double> a(7, 3, -1.0);
    for (int i = 0; i < 7; i++)
    for (int j = 0; j < 3; j++)
        ASSERT_EQ(-1.0, a(i, j));
}


TEST(Array3dTest, ValueType)
{
    ASSERT_TRUE(typeid(Array3d<double>::ValueType) == typeid(double));
}

TEST(Array3dTest, DefaultConstructor)
{
    Array3d<double> a;
    ASSERT_TRUE(Array3d<double>::IndexType(0, 0, 0) == a.getSize());
}

TEST(Array3dTest, Constructor)
{
    Array3d<double> a(4, 5, 2);
    ASSERT_TRUE(Array3d<double>::IndexType(4, 5, 2) == a.getSize());
    for (int i = 0; i < a.getSize().x; i++)
    for (int j = 0; j < a.getSize().y; j++)
    for (int k = 0; k < a.getSize().z; k++)
        ASSERT_EQ(0, a(i, j, k));
    Array3d<int> b(Array3d<double>::IndexType(5, 2, 4), 10);
    ASSERT_TRUE(Array3d<double>::IndexType(5, 2, 4) == b.getSize());
    for (int i = 0; i < b.getSize().x; i++)
    for (int j = 0; j < b.getSize().y; j++)
    for (int k = 0; k < b.getSize().z; k++)
        ASSERT_EQ(10, b(i, j, k));
}

TEST(Array3dTest, GetSize)
{
    const Array3d<float> a(17, 4, 6);
    ASSERT_TRUE(Array3d<double>::IndexType(17, 4, 6) == a.getSize());
}

TEST(Array3dTest, Access)
{
    Array3d<int> a(6, 4, 1, -1);
    for (int i = 0; i < 6; i++)
    for (int j = 0; j < 4; j++)
    for (int k = 0; k < 1; k++)
        a(i, j, k) = i + j + k;
    for (int i = 0; i < 6; i++)
    for (int j = 0; j < 4; j++)
    for (int k = 0; k < 1; k++)
        ASSERT_EQ(i + j + k, a(i, j, k));
}

TEST(Array3dTest, AccessConst)
{
    const Array3d<double> a(9, 5, 7, 3.5);
    for (int i = 0; i < 9; i++)
    for (int j = 0; j < 5; j++)
    for (int k = 0; k < 7; k++)
        ASSERT_EQ(3.5, a(i, j, k));
}
