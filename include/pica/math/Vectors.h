#ifndef PICA_VECTORS_H
#define PICA_VECTORS_H


#include "pica/math/Dimension.h"
#include "pica/math/FP.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>


namespace pica {

// Vector of 2 components of arithmetic type T (T = int, double, etc.)
// with access by .x, .y, and [], provides basic arithmetic operations
template <typename T>
struct Vector2
{
    T x, y;

    Vector2() :
        x(0), y(0) {}

    Vector2(T _x, T _y) :
        x(_x), y(_y) {}

    template<typename U>
    Vector2(const Vector2<U>& other) :
        x(other.x), y(other.y) {}

    inline T operator[](int idx) const
    {
        return *((T*)this + idx);
    }

    inline T& operator[](int idx)
    {
        return *((T*)this + idx);
    }

    inline T volume() const
    {
        return x * y;
    }

    inline T norm() const
    {
        return sqrt(x * x + y * y);
    }

    inline T norm2() const
    {
        return x * x + y * y;
    }

};

template<typename T>
inline const Vector2<T> operator + (const Vector2<T>& v1, const Vector2<T>& v2)
{
    return Vector2<T>(v1.x + v2.x, v1.y + v2.y);
}

template<typename T>
inline Vector2<T>& operator += (Vector2<T>& v1, const Vector2<T>& v2)
{
    v1.x += v2.x;
    v1.y += v2.y;
    return v1;
}

template<typename T>
inline const Vector2<T> operator - (const Vector2<T>& v1, const Vector2<T>& v2)
{
    return Vector3<T>(v1.x - v2.x, v1.y - v2.y);
}

template<typename T>
inline Vector2<T>& operator -= (Vector2<T>& v1, const Vector2<T>& v2)
{
    v1.x -= v2.x;
    v1.y -= v2.y;
    return v1;
}

template<typename T>
inline const Vector2<T> operator * (const Vector2<T>& v1, const Vector2<T>& v2)
{
    return Vector2<T>(v1.x * v2.x, v1.y * v2.y);
}

template<typename T>
inline Vector2<T>& operator *= (Vector2<T>& v1, const Vector2<T>& v2)
{
    v1.x *= v2.x;
    v1.y *= v2.y;
    return v1;
}

template<typename T>
inline const Vector2<T> operator * (const Vector2<T>& v, T a)
{
    return Vector2<T>(v.x * a, v.y * a);
}

template<typename T>
inline const Vector2<T> operator * (T a, const Vector2<T>& v)
{
    return Vector2<T>(v.x * a, v.y * a);
}

template<typename T>
inline Vector2<T>& operator *= (Vector2<T>& v, int a)
{
    v.x *= a;
    v.y *= a;
    return v;
}

template<typename T>
inline const Vector2<T> operator / (const Vector2<T>& v1, const Vector2<T>& v2)
{
    return Vector2<T>(v1.x / v2.x, v1.y / v2.y);
}

template<typename T>
inline Vector2<T>& operator /= (Vector2<T>& v1, const Vector2<T>& v2)
{
    v1.x /= v2.x;
    v1.y /= v2.y;
    return v1;
}

template<typename T>
inline const Vector2<T> operator / (const Vector2<T>& v, T a)
{
    return Vector2<T>(v.x / a, v.y / a);
}

template<typename T>
inline Vector2<T>& operator /= (Vector2<T>& v, T a)
{
    v.x /= a;
    v.y /= a;
    return v;
}

template<typename T>
inline bool operator == (const Vector2<T>& v1, const Vector2<T>& v2)
{
    return (v1.x == v2.x) && (v1.y == v2.y);
}

template<typename T>
inline bool operator != (const Vector2<T>& v1, const Vector2<T>& v2)
{
    return (v1.x != v2.x) || (v1.y != v2.y);
}

template<typename T>
inline bool operator < (const Vector2<T>& v1, const Vector2<T>& v2)
{
    return (v1.x < v2.x) && (v1.y < v2.y);
}

template<typename T>
inline bool operator <= (const Vector2<T>& v1, const Vector2<T>& v2)
{
    return (v1.x <= v2.x) && (v1.y <= v2.y);
}

template<typename T>
inline bool operator > (const Vector2<T>& v1, const Vector2<T>& v2)
{
    return (v1.x > v2.x) && (v1.y > v2.y);
}

template<typename T>
inline bool operator >= (const Vector2<T>& v1, const Vector2<T>& v2)
{
    return (v1.x >= v2.x) && (v1.y >= v2.y);
}

template<typename T>
inline std::ostream& operator<<(std::ostream& out, const Vector2<T>& v)
{
    return out << "(" << v.x << ", " << v.y << ")";
}

template<typename T>
inline T dot(const Vector2<T>& v1, const Vector2<T>& v2)
{
    return v1.x * v2.x + v1.y * v2.y;
}


// Vector of 3 components of arithmetic type T (T = int, double, etc.)
// with access by .x, .y, .z and [], provides basic arithmetic operations
template <typename T>
struct Vector3
{
    T x, y, z;

    Vector3() :
        x(0), y(0), z(0) {}

    Vector3(T _x, T _y, T _z) :
        x(_x), y(_y), z(_z) {}

    template<typename U>
    Vector3(const Vector3<U>& other) :
        x(other.x), y(other.y), z(other.z) {}

    Vector3(Vector2<T> v) :
        x(v.x), y(v.y), z(0) {}

    inline T operator[](int idx) const
    {
        return *((T*)this + idx);
    }

    inline T& operator[](int idx)
    {
        return *((T*)this + idx);
    }

    inline T volume() const
    {
        return x * y * z;
    }

    inline T norm() const
    {
        return sqrt(x * x + y * y + z * z);
    }

    inline T norm2() const
    {
        return x * x + y * y + z * z;
    }

};

template<typename T>
inline const Vector3<T> operator + (const Vector3<T>& v1, const Vector3<T>& v2)
{
    return Vector3<T>(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

template<typename T>
inline Vector3<T>& operator += (Vector3<T>& v1, const Vector3<T>& v2)
{
    v1.x += v2.x;
    v1.y += v2.y;
    v1.z += v2.z;
    return v1;
}

template<typename T>
inline const Vector3<T> operator - (const Vector3<T>& v1, const Vector3<T>& v2)
{
    return Vector3<T>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

template<typename T>
inline Vector3<T>& operator -= (Vector3<T>& v1, const Vector3<T>& v2)
{
    v1.x -= v2.x;
    v1.y -= v2.y;
    v1.z -= v2.z;
    return v1;
}

template<typename T>
inline const Vector3<T> operator * (const Vector3<T>& v1, const Vector3<T>& v2)
{
    return Vector3<T>(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
}

template<typename T>
inline Vector3<T>& operator *= (Vector3<T>& v1, const Vector3<T>& v2)
{
    v1.x *= v2.x;
    v1.y *= v2.y;
    v1.z *= v2.z;
    return v1;
}

template<typename T>
inline const Vector3<T> operator * (const Vector3<T>& v, T a)
{
    return Vector3<T>(v.x * a, v.y * a, v.z * a);
}

template<typename T>
inline const Vector3<T> operator * (T a, const Vector3<T>& v2)
{
    return Vector3<T>(a * v2.x, a * v2.y, a * v2.z);
}

template<typename T>
inline Vector3<T>& operator *= (Vector3<T>& v, T a)
{
    v.x *= a;
    v.y *= a;
    v.z *= a;
    return v;
}

template<typename T>
inline const Vector3<T> operator / (const Vector3<T>& v1, const Vector3<T>& v2)
{
    return Vector3<T>(v1.x / v2.x, v1.y / v2.y, v1.z / v2.z);
}

template<typename T>
inline Vector3<T>& operator /= (Vector3<T>& v1, const Vector3<T>& v2)
{
    v1.x /= v2.x;
    v1.y /= v2.y;
    v1.z /= v2.z;
    return v1;
}

template<typename T>
inline const Vector3<T> operator / (const Vector3<T>& v, T a)
{
    return Vector3<T>(v.x / a, v.y / a, v.z / a);
}

template<typename T>
inline Vector3<T>& operator /= (Vector3<T>& v, T a)
{
    v.x /= a;
    v.y /= a;
    v.z /= a;
    return v;
}

template<typename T>
inline bool operator == (const Vector3<T>& v1, const Vector3<T>& v2)
{
    return (v1.x == v2.x) && (v1.y == v2.y) && (v1.z == v2.z);
}

template<typename T>
inline bool operator != (const Vector3<T>& v1, const Vector3<T>& v2)
{
    return (v1.x != v2.x) || (v1.y != v2.y) || (v1.z != v2.z);
}

template<typename T>
inline bool operator < (const Vector3<T>& v1, const Vector3<T>& v2)
{
    return (v1.x < v2.x) && (v1.y < v2.y) && (v1.z < v2.z);
}

template<typename T>
inline bool operator <= (const Vector3<T>& v1, const Vector3<T>& v2)
{
    return (v1.x <= v2.x) && (v1.y <= v2.y) && (v1.z <= v2.z);
}

template<typename T>
inline bool operator > (const Vector3<T>& v1, const Vector3<T>& v2)
{
    return (v1.x > v2.x) && (v1.y > v2.y) && (v1.z > v2.z);
}

template<typename T>
inline bool operator >= (const Vector3<T>& v1, const Vector3<T>& v2)
{
    return (v1.x >= v2.x) && (v1.y >= v2.y) && (v1.z >= v2.z);
}

template<typename T>
inline std::ostream& operator<<(std::ostream& out, const Vector3<T>& v)
{
    return out << "(" << v.x << ", " << v.y << ", " << v.z << ")";
}

template<typename T>
inline T dot(const Vector3<T>& v1, const Vector3<T>& v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

template<typename T>
inline const Vector3<T> cross(const Vector3<T>& v1, const Vector3<T>& v2)
{
    return Vector3<T>(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z,
        v1.x * v2.y - v1.y * v2.x);
}


typedef Vector3<int> Int3;
typedef Vector3<FP> FP3;

inline const Int3 truncate(const FP3& v)
{
    return Int3((int)v.x, (int)v.y, (int)v.z);
}

inline Int3 remainder(const Int3& v1, const Int3& v2)
{
    return Int3(v1.x % v2.x, v1.y % v2.y, v1.z % v2.z);
}

inline FP dist(const FP3& v1, const FP3& v2)
{
    return (v1 - v2).norm();
}


template<Dimension dimension>
struct VectorFPTraits {
};

template<>
struct VectorFPTraits<One> {
    typedef FP Type;
};

template<>
struct VectorFPTraits<Two> {
    typedef Vector2<FP> Type;
};

template<>
struct VectorFPTraits<Three> {
    typedef Vector3<FP> Type;
};


} // namespace pica


#endif
