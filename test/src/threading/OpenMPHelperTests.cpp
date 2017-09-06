#include "TestingUtility.h"

#include "pica/threading/OpenMPHelper.h"

#include <vector>

using namespace pica;


/// For an unknown yet reason this test does not compile on mvs100k
/*
TEST(UtilsOpenMPHelperTest, GetNumThreads) {
    std::vector<int> numThreadsVector(getNumThreads());
    int numThreads = (int)numThreadsVector.size();
    #pragma omp parallel for
    for (int i = 0; i < numThreads; i++)
        numThreadsVector[i] = omp_get_num_threads();
    for (int i = 0; i < numThreads; i++)
        ASSERT_EQ(numThreads, numThreadsVector[i]);
}*/

TEST(OpenMPHelperTest, UseOpenMP) {
    bool enabled = useOpenMP();
    if (getNumThreads() > 1)
        ASSERT_TRUE(enabled);
}

TEST(OpenMPHelperTest, OpenMPDisabled) {
    if (!useOpenMP()) {
        ASSERT_EQ(1, omp_get_max_threads());
        ASSERT_EQ(0, omp_get_thread_num());
        ASSERT_EQ(1, omp_get_num_threads());
    }
}

TEST(OpenMPHelperTest, Locks) {
    // test that lock functions do not crash
    omp_lock_t lock;
    omp_init_lock(&lock);
    omp_set_lock(&lock);
    omp_unset_lock(&lock);
    omp_destroy_lock(&lock);
}
