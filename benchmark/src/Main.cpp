#include "pica/utility/Utility.h"
#include "pica/threading/OpenMPHelper.h"

#include "Parameters.h"

#include <iostream>
#include <string>

using namespace pica;
using namespace std;


int realMain(int argc, char* argv[]);

int main(int argc, char* argv[])
{
    try {
        return realMain(argc, argv);
    }
    catch (...) {
        cout << "ERROR: Unhandled exception\n";
        return EXIT_FAILURE;
    }
}

int realMain(int argc, char* argv[])
{
    string message = "pica benchmark, ";
    if (useOpenMP())
        message += toString(getNumThreads()) + " OpenMP threads.";
    else
        message += "OpenMP disabled.";
    cout << message << "\n\n";

    Parameters parameters;

    return 0;
}
