#include "pica/Utility.h"

#include <iostream>
#include <limits>
#include <sstream>


namespace pica {

std::string toString(int i)
{
    std::stringstream ss;
    ss << i;
    return ss.str();
}

std::string toString(double x)
{
    std::stringstream ss;
    ss << x;
    return ss.str();
}

bool isNumberFinite(double x)
{
    return x <= std::numeric_limits<double>::max()
        && x >= -std::numeric_limits<double>::max();
}

} // namespace pica
