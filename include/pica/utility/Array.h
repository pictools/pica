#ifndef PICA_ARRAY_H
#define PICA_ARRAY_H


#include "pica/utility/Assert.h"
#include "pica/math/Dimension.h"
#include "pica/math/Vectors.h"

#include <vector>


namespace pica {


// A simple class for 1d array (vector) of type T
template<typename T>
class Array1d {
public:
    typedef T ValueType;
    typedef Vector1<int> IndexType;

    // Create an empty vector
    Array1d() {}

    // Create an array of the given size with the given value of elements
    Array1d(IndexType size, ValueType value = 0) :
        data(size.x, value) {}

    // Get size of the vector
    IndexType getSize() const { return IndexType(static_cast<int>(data.size())); }

    // Access the value by index without checking the index
    ValueType& operator()(int i)
    {
        PICA_ASSERT((i >= 0) && (i < getSize().x));
        return data[i];
    }

    const ValueType& operator()(int i) const
    {
        PICA_ASSERT((i >= 0) && (i < getSize().x));
        return data[i];
    }

    ValueType& operator()(IndexType index)
    {
        PICA_ASSERT((index.x >= 0) && (index.x < getSize().x));
        return data[index.x];
    }

    const ValueType& operator()(IndexType index) const
    {
        PICA_ASSERT((index.x >= 0) && (index.x < getSize().x));
        return data[index.x];
    }

private:
    std::vector<T> data;
};


// A simple class for 2d array (vector) of type T
template<typename T>
class Array2d {
public:
    typedef T ValueType;
    typedef Vector2<int> IndexType;

    // Create an empty matrix
    Array2d() :
        size(0, 0) {}

    // Create a matrix of the given size with the given value of elements
    Array2d(int nRows, int nCols, ValueType value = 0) :
        size(IndexType(nRows, nCols)), data(size.volume(), value) {}

    Array2d(IndexType size, ValueType value = 0) :
        size(size), data(size.volume(), value) {}

    // Get number of rows, columns and total number of elements
    IndexType getSize() const { return size; }
    int getNumRows() const { return size.x; }
    int getNumCols() const { return size.y; }
    int getNumElements() const { return (int)data.size(); }

    // Access the value by index without checking the index
    ValueType& operator()(int i, int j)
    {
        PICA_ASSERT((i >= 0) && (i < size.x) && (j >= 0) && (j < size.y));
        return data[i * size.y + j];
    }

    const ValueType& operator()(int i, int j) const
    {
        PICA_ASSERT((i >= 0) && (i < size.x) && (j >= 0) && (j < size.y));
        return data[i * size.y + j];
    }

    ValueType& operator()(IndexType index)
    {
        PICA_ASSERT((index.x >= 0) && (index.x < size.x) && (index.y >= 0) && (index.y < size.y));
        return data[index.x * size.y + index.y];
    }

    const ValueType& operator()(IndexType index) const
    {
        PICA_ASSERT((index.x >= 0) && (index.x < size.x) && (index.y >= 0) && (index.y < size.y));
        return data[index.x * size.y + index.y];
    }

private:
    IndexType size;
    std::vector<T> data;
};


// A simple class for 3d array (vector) of type T
template<typename T>
class Array3d  {
public:
    typedef T ValueType;
    typedef Vector3<int> IndexType;

    // Create an empty 3d array
    Array3d() :
        size(0, 0, 0) {}

    // Create a matrix of the given size with the given value of elements
    Array3d(int n1, int n2, int n3, T value = 0) :
        size(IndexType(n1, n2, n3)), data(size.x * size.y * size.z, value) {}

    Array3d(IndexType size, T value = 0) :
        size(size), data(size.x * size.y * size.z, value) {}

    // Get size and total number of elements
    IndexType getSize() const { return size; }
    int getNumElements() const { return data.size(); }

    // Access the value by index without checking the index
    ValueType& operator()(int i, int j, int k)
    {
        PICA_ASSERT((i >= 0) && (i < size.x) && (j >= 0) && (j < size.y) && (k >= 0) && (k < size.z));
        return data[(i * size.y + j) * size.z + k];
    }

    const ValueType& operator()(int i, int j, int k) const
    {
        PICA_ASSERT((i >= 0) && (i < size.x) && (j >= 0) && (j < size.y) && (k >= 0) && (k < size.z));
        return data[(i * size.y + j) * size.z + k];
    }

    ValueType& operator()(IndexType index)
    {
        PICA_ASSERT((index.x >= 0) && (index.x < size.x) && (index.y >= 0) && (index.y < size.y) && (index.z >= 0) && (index.z < size.z));
        return data[(index.x * size.y + index.y) * size.z + index.z];
    }

    const ValueType& operator()(IndexType index) const
    {
        PICA_ASSERT((index.x >= 0) && (index.x < size.x) && (index.y >= 0) && (index.y < size.y) && (index.z >= 0) && (index.z < size.z));
        return data[(index.x * size.y + index.y) * size.z + index.z];
    }

private:
    IndexType size;
    std::vector<T> data;
};


template<Dimension dimension, typename T>
struct ArrayTypeHelper {
};

template<typename T>
struct ArrayTypeHelper<One, T> {
    typedef Array1d<T> Type;
};

template<typename T>
struct ArrayTypeHelper<Two, T> {
    typedef Array2d<T> Type;
};

template<typename T>
struct ArrayTypeHelper<Three, T> {
    typedef Array3d<T> Type;
};


template<int dimension, typename T>
struct ArrayIntTypeHelper {
};

template<typename T>
struct ArrayIntTypeHelper<1, T> {
    typedef Array1d<T> Type;
};

template<typename T>
struct ArrayIntTypeHelper<2, T> {
    typedef Array2d<T> Type;
};

template<typename T>
struct ArrayIntTypeHelper<3, T> {
    typedef Array3d<T> Type;
};


} // namespace pica


#endif
