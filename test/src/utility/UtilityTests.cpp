#include "TestingUtility.h"

#include "pica/utility/Utility.h"

#include <limits>
#include <string>
#include <vector>
using std::string;
using std::vector;

using namespace pica;


TEST(UtilityTest, Ptr) {
    vector<int> vInt;
    vInt.push_back(3);
    vInt.push_back(-6);
    vInt.push_back(4);
    vInt.push_back(0);
    int* ptrInt = ptr(vInt);
    for (int i = 0; i < vInt.size(); i++)
        ASSERT_EQ(vInt[i], ptrInt[i]);
    vector<double> emptyVector;
    ASSERT_EQ(NULL, ptr(emptyVector));
}

TEST(UtilityTest, PtrConst) {
    const vector<double> vDouble(5, 1.3);
    const double* ptrDouble = ptr(vDouble);
    for (int i = 0; i < vDouble.size(); i++)
        ASSERT_EQ(vDouble[i], ptrDouble[i]);
    const vector<char> emptyVector;
    ASSERT_EQ(NULL, ptr(emptyVector));
}

TEST(UtilityTest, IsNumberFinite) {
    ASSERT_TRUE(isNumberFinite(0.0));
    ASSERT_TRUE(isNumberFinite(1.31));
    ASSERT_TRUE(isNumberFinite(-1.12e-5));
    ASSERT_TRUE(isNumberFinite(std::numeric_limits<double>::denorm_min()));
    ASSERT_TRUE(isNumberFinite(-std::numeric_limits<double>::denorm_min()));
    ASSERT_TRUE(isNumberFinite(std::numeric_limits<double>::max()));
    ASSERT_TRUE(isNumberFinite(-std::numeric_limits<double>::max()));
    ASSERT_FALSE(isNumberFinite(std::numeric_limits<double>::infinity()));
    ASSERT_FALSE(isNumberFinite(-std::numeric_limits<double>::infinity()));
    ASSERT_FALSE(isNumberFinite(std::numeric_limits<double>::quiet_NaN()));
    ASSERT_FALSE(isNumberFinite(std::numeric_limits<double>::signaling_NaN()));
}


TEST(UtilityTest, IntToString)
{
    ASSERT_EQ("0", toString(0));
    ASSERT_EQ("5", toString(5));
    ASSERT_EQ("23621", toString(23621));
    ASSERT_EQ("-72", toString(-72));
    ASSERT_EQ("-563534", toString(-563534));
}

TEST(UtilityTest, DoubleToString)
{
    ASSERT_EQ("6.5", toString(6.5));
    ASSERT_EQ("-13.25", toString(-13.25));
}

TEST(UtilityTest, BoolToString)
{
    ASSERT_EQ("1", toString(true));
    ASSERT_EQ("0", toString(false));
}

TEST(UtilityTest, Int3ToString)
{
    ASSERT_EQ("(2, -1, 0)", toString(Int3(2, -1, 0)));
    ASSERT_EQ("(-5  xyz -4  xyz 6)", toString(Int3(-5, -4, 6), "  xyz "));
    ASSERT_EQ("[ 3-5-2)", toString(Int3(3, 5, 2), "-", "[ "));
    ASSERT_EQ("(15,4,-3]", toString(Int3(15, 4, -3), ",", "(", "]"));
}

TEST(UtilityTest, FP3ToString)
{
    ASSERT_EQ("(2, -0.25, 1.5)", toString(FP3(2, -0.25, 1.5)));
    ASSERT_EQ("(2  xyz -0.25  xyz 1.5)", toString(FP3(2, -0.25, 1.5), "  xyz "));
    ASSERT_EQ("[ 2-0.25-1.5)", toString(FP3(2, 0.25, 1.5), "-", "[ "));
    ASSERT_EQ("(2,-0.25,1.5]", toString(FP3(2, -0.25, 1.5), ",", "(", "]"));
}

TEST(IndexingTest, GetLinearIndex)
{
    ASSERT_EQ(6 * 2 * 2 + 2 * 3 + 1, getLinearIndex(Int3(4, 6, 2), Int3(2, 3, 1)));
    ASSERT_EQ(5 * 8 * 7 + 8 * 4 + 2, getLinearIndex(Int3(12, 5, 8), Int3(7, 4, 2)));
}

TEST(IndexingTest, getReorderedLinearIndex)
{
    ASSERT_EQ(6 * 4 * 1 + 4 * 3 + 2, getReorderedLinearIndex(Int3(4, 6, 2), Int3(2, 3, 1)));
    ASSERT_EQ(5 * 12 * 2 + 12 * 4 + 7, getReorderedLinearIndex(Int3(12, 5, 8), Int3(7, 4, 2)));
}


class IndexIntervalTest : public testing::Test {
protected:
    virtual void SetUp()
    { 
        begin = Int3(1, 2, 4);
        end = Int3(5, 12, 7);
        step = Int3(2, 3, 1);
    }
    Int3 begin, end, step;
};

TEST_F(IndexIntervalTest, ValuesSetInConstructor)
{
    IndexInterval<Int3> interval(begin, end, step);
    ASSERT_EQ_INT3(begin, interval.begin);
    ASSERT_EQ_INT3(end, interval.end);
    ASSERT_EQ_INT3(step, interval.step);
}

TEST_F(IndexIntervalTest, DefaultStepValue)
{
    IndexInterval<Int3> interval(begin, end);
    ASSERT_EQ_INT3(begin, interval.begin);
    ASSERT_EQ_INT3(end, interval.end);
    ASSERT_EQ_INT3(Int3(1, 1, 1), interval.step);
}

TEST_F(IndexIntervalTest, Volume)
{
    IndexInterval<Int3> interval(begin, end, step);
    // 2 * 4 * 3 = 24 indexes in the index space
    ASSERT_EQ(24, interval.volume());
}

TEST_F(IndexIntervalTest, IncrementIndex)
{
    IndexInterval<Int3> interval(begin, end, step);
    Int3 idx = begin;
    for (int i = begin.x; i < end.x; i += step.x)
    for (int j = begin.y; j < end.y; j += step.y)
    for (int k = begin.z; k < end.z; k += step.z) {
        ASSERT_EQ_INT3(Int3(i, j, k), idx);
        interval.incrementIndex(idx);
    }
}

TEST_F(IndexIntervalTest, ToInt)
{
    IndexInterval<Int3> interval(begin, end, step);
    int idx = 0;
    for (int i = begin.x; i < end.x; i += step.x)
    for (int j = begin.y; j < end.y; j += step.y)
    for (int k = begin.z; k < end.z; k += step.z) {
        ASSERT_EQ(idx, interval.toInt(Int3(i, j, k)));
        idx++;
    }
}

TEST_F(IndexIntervalTest, ToIndex3d)
{
    IndexInterval<Int3> interval(begin, end, step);
    int idx = 0;
    for (int i = begin.x; i < end.x; i += step.x)
    for (int j = begin.y; j < end.y; j += step.y)
    for (int k = begin.z; k < end.z; k += step.z) {
        ASSERT_EQ_INT3(Int3(i, j, k), interval.toIndex3d(idx));
        idx++;
    }
}