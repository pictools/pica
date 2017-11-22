#ifndef PICA_UTILITY_H
#define PICA_UTILITY_H


#include <limits>
#include <string>
#include <sstream>
#include <vector>


namespace pica {

/* Return a pointer to the first element of a vector. Note that the elements
   are required by the standard to be stored contiguously. */
template<typename T>
T * ptr(std::vector<T> & v)
{
    if (v.empty()) return NULL;
    return &v.front();
}

template<typename T>
const T * ptr(const std::vector<T> & v)
{
    if (v.empty()) return NULL;
    return &v.front();
}

// Check if a number is finite (e.g. not an infinity or NaN). Poor man's
// isfinite() from C99.
inline bool isNumberFinite(double x)
{
    return x <= std::numeric_limits<double>::max()
        && x >= -std::numeric_limits<double>::max();
}

// Convert to string value of type T which can be written to ostream 
template<typename T>
inline std::string toString(const T& value)
{
    std::ostringstream stream;
    stream << value;
    return stream.str();
}


// Version for vectors of 3 components with .x, .y, .z access
template<class Vector3>
inline std::string toString(const Vector3& v, const std::string& delimeter,
    const std::string& prefix = "(", const std::string& postfix = ")")
{
    std::ostringstream stream;
    stream << prefix << v.x << delimeter << v.y << delimeter << v.z << postfix;
    return stream.str();
}


template<typename Index3d>
inline int getLinearIndex(const Index3d& size, const Index3d& index)
{
    return size.y * size.z * index.x + size.z * index.y + index.z;
}

template<typename Index3d>
inline int getReorderedLinearIndex(const Index3d& size, const Index3d& index)
{
    return size.x * size.y * index.z + size.x * index.y + index.x;
}


/* Index space for half-closed 3d interval [begin, end) with positive step.
Uses x > y > z order. Supports increment and conversion between 3d indexes
in the interval and linear 0-based indexes. */
template<typename Index3d>
struct IndexInterval {
    const Index3d begin, end, step, size, shiftedBegin;

    IndexInterval(const Index3d& _begin, const Index3d& _end,
                  const Index3d& _step = Index3d(1, 1, 1)):
        begin(_begin), end(_end), step(_step),
        size((_end - _begin + _step - Index3d(1, 1, 1)) / _step),
        shiftedBegin(_begin - _step + Index3d(1, 1, 1)) {
    }

    int volume() const {
        return size.volume();
    }

    void incrementIndex(Index3d& idx) const {
        for (int d = 2; d >= 0; --d)
            if (idx[d] < end[d] - step[d]) {
                idx[d] += step[d];
                return;
            }
            else
                idx[d] = begin[d];
        // if an index can not be incremented, set to end
        idx = end;
    }

    int toInt(const Index3d& idx) const
    {
        Index3d diff = (idx - shiftedBegin) / step;
        return diff.z + size.z * (diff.y + size.y * diff.x);
    }

    const Index3d toIndex3d(int idx) const
    {
        return begin + Index3d(idx / (size.z * size.y),
            (idx / size.z) % size.y, idx % size.z) * step;
    }

};

} // namespace pica


#endif
