#include "Output.h"
#include "Parameters.h"
#include "Parser.h"
#include "PerformanceTracker.h"

#include <iostream>
#include <map>
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
        cout << "ERROR: Unhandled exception, benchmark terminated\n";
        return EXIT_FAILURE;
    }
}

int realMain(int argc, char* argv[])
{
    printHeader();
    Parameters parameters = readParameters(argc, argv);
    printParameters(parameters);
    PerformanceTracker tracker;
    printPerformance(tracker);
    return 0;
}
