#include "TestingUtility.h"

#include "pica/grid/YeeGrid.h"

#include <memory>

using namespace pica;


class YeeGridTest1d : public ::testing::Test {
public:
    virtual void SetUp()
    {
        origin = -1.5;
        step = 0.5;
        size = 10;
        grid.reset(new GridType(origin, step, size));
    }
    typedef YeeGrid<One, double> GridType;
    typedef GridType::PositionType PositionType;
    typedef GridType::IndexType IndexType;
    std::auto_ptr<GridType> grid;
    PositionType origin, step;
    IndexType size;
};

TEST_F(YeeGridTest1d, GetOrigin)
{
    EXPECT_TRUE(origin == grid->getOrigin());
}

TEST_F(YeeGridTest1d, GetStep)
{
    EXPECT_TRUE(step == grid->getStep());
}

TEST_F(YeeGridTest1d, GetSize)
{
    EXPECT_TRUE(size == grid->getSize());
}

TEST_F(YeeGridTest1d, GetCellIndex)
{
    for (int i = 0; i < grid->getSize().x; i++) {
        PositionType position = grid->getOrigin() + grid->getStep() * (i + 0.5);
        EXPECT_TRUE(IndexType(i) == grid->getCellIndex(position));
    }
}


class YeeGridTest2d : public ::testing::Test {
public:
    virtual void SetUp()
    {
        origin = PositionType(3.5, -8.6);
        step = PositionType(0.1, 0.2);
        size = IndexType(5, 6);
        grid.reset(new GridType(origin, step, size));
    }
    typedef YeeGrid<Two, double> GridType;
    typedef GridType::PositionType PositionType;
    typedef GridType::IndexType IndexType;
    std::auto_ptr<GridType> grid;
    PositionType origin, step;
    IndexType size;
};

TEST_F(YeeGridTest2d, GetOrigin)
{
    EXPECT_TRUE(origin == grid->getOrigin());
}

TEST_F(YeeGridTest2d, GetStep)
{
    EXPECT_TRUE(step == grid->getStep());
}

TEST_F(YeeGridTest2d, GetSize)
{
    EXPECT_TRUE(size == grid->getSize());
}

TEST_F(YeeGridTest2d, GetCellIndex)
{
    for (int i = 0; i < grid->getSize().x; i++)
    for (int j = 0; j < grid->getSize().y; j++)
    {
        PositionType position = grid->getOrigin() + grid->getStep() * PositionType(i + 0.5, j + 0.5);
        EXPECT_TRUE(IndexType(i, j) == grid->getCellIndex(position));
    }
}


class YeeGridTest3d : public ::testing::Test {
public:
    virtual void SetUp()
    {
        origin = PositionType(-1.5, 3.4, 0.31);
        step = PositionType(0.1, 0.2, 0.15);
        size = IndexType(12, 3, 4);
        grid.reset(new GridType(origin, step, size));
    }
    typedef YeeGrid<Three, double> GridType;
    typedef GridType::PositionType PositionType;
    typedef GridType::IndexType IndexType;
    std::auto_ptr<GridType> grid;
    PositionType origin, step;
    IndexType size;
};

TEST_F(YeeGridTest3d, GetOrigin)
{
    EXPECT_TRUE(origin == grid->getOrigin());
}

TEST_F(YeeGridTest3d, GetStep)
{
    EXPECT_TRUE(step == grid->getStep());
}

TEST_F(YeeGridTest3d, GetSize)
{
    EXPECT_TRUE(size == grid->getSize());
}

TEST_F(YeeGridTest3d, GetCellIndex)
{
    for (int i = 0; i < grid->getSize().x; i++)
    for (int j = 0; j < grid->getSize().y; j++)
    for (int k = 0; k < grid->getSize().z; k++) {
        PositionType position = grid->getOrigin() + grid->getStep() * PositionType(i + 0.5, j + 0.5, k + 0.5);
        EXPECT_TRUE(IndexType(i, j, k) == grid->getCellIndex(position));
    }
}
