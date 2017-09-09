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
    return Vector2<T>(v1.x - v2.x, v1.y - v2.y);
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


template<Dimension dimension, typename T>
struct VectorTypeHelper {
};

template<typename T>
struct VectorTypeHelper<One, T> {
    typedef T Type;
};

template<typename T>
struct VectorTypeHelper<Two, T> {
    typedef Vector2<T> Type;
};

template<typename T>
struct VectorTypeHelper<Three, T> {
    typedef Vector3<T> Type;
};

typedef Vector3<int> Int3;
typedef Vector3<FP> FP3;

template<typename Real>
inline const int truncate(const Real& v)
{
    return (int)v;
}

template<typename Real>
inline const Vector2<int> truncate(const Vector2<Real>& v)
{
    return Vector2<int>(static_cast<int>(v.x), static_cast<int>(v.y));
}

template<typename Real>
inline const Vector3<int> truncate(const Vector3<Real>& v)
{
    return Vector3<int>(static_cast<int>(v.x), static_cast<int>(v.y), static_cast<int>(v.z));
}

template<typename Real>
inline Real inverse(const Real a)
{
    return static_cast<Real>(1) / a;
}

template<typename Real>
inline const Vector2<Real> inverse(const Vector2<Real>& v)
{
    return Vector2<Real>(static_cast<Real>(1) / v.x, static_cast<Real>(1) / v.y);
}

template<typename Real>
inline const Vector3<Real> inverse(const Vector3<Real>& v)
{
    return Vector3<Real>(static_cast<Real>(1) / v.x, static_cast<Real>(1) / v.y, static_cast<Real>(1) / v.z);
}

template<Dimension dimension, typename T>
struct OnesHelper {
    static typename VectorTypeHelper<dimension, T>::Type get();
};

template<typename T>
struct OnesHelper<One, T> {
    static typename VectorTypeHelper<One, T>::Type get() { return static_cast<T>(1); }
};

template<typename T>
struct OnesHelper<Two, T> {
    static typename VectorTypeHelper<Two, T>::Type get() { return VectorTypeHelper<Two, T>::Type(static_cast<T>(1), static_cast<T>(1)); }
};

template<typename T>
struct OnesHelper<Three, T> {
    static typename VectorTypeHelper<Three, T>::Type get() { return VectorTypeHelper<Three, T>::Type(static_cast<T>(1), static_cast<T>(1), static_cast<T>(1)); }
};

template<Dimension dimension, typename T>
inline typename VectorTypeHelper<dimension, T>::Type ones()
{
    return OnesHelper<dimension, T>::get();
}


template<typename T>
struct ScalarType {
    typedef T Type;
};

template<typename T>
struct ScalarType<Vector2<T> > {
    typedef T Type;
};

template<typename T>
struct ScalarType<Vector3<T> > {
    typedef T Type;
};

inline Int3 remainder(const Int3& v1, const Int3& v2)
{
    return Int3(v1.x % v2.x, v1.y % v2.y, v1.z % v2.z);
}

inline FP dist(const FP3& v1, const FP3& v2)
{
    return (v1 - v2).norm();
}





} // namespace pica


#endif
