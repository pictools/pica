#include "TestingUtility.h"

#include "pica/threading/OpenMPHelper.h"
#include "pica/utility/Utility.h"

#include <cstdio>
#include <iostream>
#include <string>


int main(int argc, char **argv)
{
    /* When running under MPI, stderr becomes redirected and buffered.
       We need it to be unbuffered, so that if we abort, any messages get
       through first. */
    std::setvbuf(stderr, NULL, _IONBF, 0);
    std::ostream& out = std::cout;
    std::ostream& err = std::cerr;

    out << "Running tests built on " << __DATE__ << " " << __TIME__ << "\n";
    #ifdef __MIC
        out << "Using Xeon Phi, ";
    #else
        out << "Using CPU, ";
    #endif
    if (pica::useOpenMP())
        out << pica::getNumThreads() << " OpenMP threads.\n\n";
    else
        out << " OpenMP disabled.\n\n";
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

    return result;
}
