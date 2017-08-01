#ifndef PICA_VECTORS_H
#define PICA_VECTORS_H


#include "pica/FP.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>


namespace pica {

// Vector of 3 ints. Can be accessed by either .x, .y, .z or [0], [1], [2].
struct Int3
{
    int x, y, z;

    Int3():
        x(0), y(0), z(0) {}

    Int3(int _x, int _y, int _z):
        x(_x), y(_y), z(_z) {}

    inline int operator[](int idx) const
    { return *((int*)this + idx); }

    inline int& operator[](int idx)
    { return *((int*)this + idx); }

    inline int volume() const
    { return x * y * z; }
};

inline const Int3 operator + (const Int3& v1, const Int3& v2)
{
    return Int3(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

inline Int3& operator += (Int3& v1, const Int3& v2)
{
    v1.x += v2.x;
    v1.y += v2.y;
    v1.z += v2.z;
    return v1;
}

inline const Int3 operator - (const Int3& v1, const Int3& v2)
{
    return Int3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

inline Int3& operator -= (Int3& v1, const Int3& v2)
{
    v1.x -= v2.x;
    v1.y -= v2.y;
    v1.z -= v2.z;
    return v1;
}

inline const Int3 operator * (const Int3& v1, const Int3& v2)
{
    return Int3(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
}

inline Int3& operator *= (Int3& v1, const Int3& v2)
{
    v1.x *= v2.x;
    v1.y *= v2.y;
    v1.z *= v2.z;
    return v1;
}

inline const Int3 operator * (const Int3& v, int a)
{
    return Int3(v.x * a, v.y * a, v.z * a);
}

inline Int3& operator *= (Int3& v, int a)
{
    v.x *= a;
    v.y *= a;
    v.z *= a;
    return v;
}

inline const Int3 operator / (const Int3& v1, const Int3& v2)
{
    return Int3(v1.x / v2.x, v1.y / v2.y, v1.z / v2.z);
}

inline Int3& operator /= (Int3& v1, const Int3& v2)
{
    v1.x /= v2.x;
    v1.y /= v2.y;
    v1.z /= v2.z;
    return v1;
}

inline const Int3 operator / (const Int3& v, int a)
{
    return Int3(v.x / a, v.y / a, v.z / a);
}

inline Int3& operator /= (Int3& v, int a)
{
    v.x /= a;
    v.y /= a;
    v.z /= a;
    return v;
}

inline bool operator == (const Int3& v1, const Int3& v2)
{
    return (v1.x == v2.x) && (v1.y == v2.y) && (v1.z == v2.z);
}

inline bool operator != (const Int3& v1, const Int3& v2)
{
    return (v1.x != v2.x) || (v1.y != v2.y) || (v1.z != v2.z);
}

inline bool operator < (const Int3& v1, const Int3& v2)
{
    return (v1.x < v2.x) && (v1.y < v2.y) && (v1.z < v2.z);
}

inline bool operator <= (const Int3& v1, const Int3& v2)
{
    return (v1.x <= v2.x) && (v1.y <= v2.y) && (v1.z <= v2.z);
}

inline bool operator > (const Int3& v1, const Int3& v2)
{
    return (v1.x > v2.x) && (v1.y > v2.y) && (v1.z > v2.z);
}

inline bool operator >= (const Int3& v1, const Int3& v2)
{
    return (v1.x >= v2.x) && (v1.y >= v2.y) && (v1.z >= v2.z);
}

inline Int3 remainder(const Int3& v1, const Int3& v2)
{
    return Int3(v1.x % v2.x, v1.y % v2.y, v1.z % v2.z);
}

inline std::ostream& operator<<(std::ostream& out, const Int3& v)
{
    return out << "(" << v.x << ", " << v.y << ", " << v.z << ")";
}


// Vector of 3 FPs. Can be accessed by either .x, .y, .z or [0], [1], [2].
struct FP3
{
    FP x, y, z;

    FP3():
        x(0), y(0), z(0) {}

    FP3(FP _x, FP _y, FP _z):
        x(_x), y(_y), z(_z) {}

    explicit FP3(const Int3& v):
        x(v.x), y(v.y), z(v.z) {}

    inline FP operator[](int idx) const
    { return *((FP*)this + idx); }

    inline FP& operator[](int idx)
    { return *((FP*)this + idx); }

    inline FP norm() const
    { return sqrt(x * x + y * y + z * z); }

    inline FP norm2() const
    { return x * x + y * y + z * z; }

    inline FP volume() const
    { return x * y * z; }
};

inline const FP3 operator + (const FP3& v1, const FP3& v2)
{
    return FP3(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

inline FP3& operator += (FP3& v1, const FP3& v2)
{
    v1.x += v2.x;
    v1.y += v2.y;
    v1.z += v2.z;
    return v1;
}

inline const FP3 operator - (const FP3& v1, const FP3& v2)
{
   return FP3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

inline FP3& operator -= (FP3& v1, const FP3& v2)
{
    v1.x -= v2.x;
    v1.y -= v2.y;
    v1.z -= v2.z;
    return v1;
}

inline const FP3 operator * (const FP3& v1, const FP3& v2)
{
    return FP3(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
}

inline FP3& operator *= (FP3& v1, const FP3& v2)
{
    v1.x *= v2.x;
    v1.y *= v2.y;
    v1.z *= v2.z;
    return v1;
}

inline const FP3 operator * (const FP3& v1, const Int3& v2)
{
   return FP3(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
}

inline FP3& operator *= (FP3& v1, const Int3& v2)
{
    v1.x *= v2.x;
    v1.y *= v2.y;
    v1.z *= v2.z;
    return v1;
}

inline const FP3 operator * (const FP3& v, FP a)
{
    return FP3(v.x * a, v.y * a, v.z * a);
}

inline FP3& operator *= (FP3& v, FP a)
{
    v.x *= a;
    v.y *= a;
    v.z *= a;
    return v;
}

inline const FP3 operator * (FP a, const FP3& v1)
{
    return FP3(v1.x * a, v1.y * a, v1.z * a);
}

inline const FP3 operator / (const FP3 & v1, const FP3 & v2)
{
    return FP3(v1.x / v2.x, v1.y / v2.y, v1.z / v2.z);
}

inline FP3& operator /= (FP3& v1, const FP3& v2)
{
    v1.x /= v2.x;
    v1.y /= v2.y;
    v1.z /= v2.z;
    return v1;
}

inline const FP3 operator / (const FP3 & v1, const Int3 & v2)
{
    return FP3(v1.x / v2.x, v1.y / v2.y, v1.z / v2.z);
}

inline FP3& operator /= (FP3& v1, const Int3& v2)
{
    v1.x /= v2.x;
    v1.y /= v2.y;
    v1.z /= v2.z;
    return v1;
}

inline const FP3 operator / (const FP3 & v, FP a)
{
    return FP3(v.x / a, v.y / a, v.z / a);
}

inline FP3& operator /= (FP3& v, FP a)
{
    v.x /= a;
    v.y /= a;
    v.z /= a;
    return v;
}

inline bool operator == (const FP3 & v1, const FP3 & v2)
{
    return (v1.x == v2.x) && (v1.y == v2.y) && (v1.z == v2.z);
}

inline bool operator != (const FP3 & v1, const FP3 & v2)
{
    return (v1.x != v2.x) || (v1.y != v2.y) || (v1.z != v2.z);
}

inline bool operator < (const FP3 & v1, const FP3 & v2)
{
    return (v1.x < v2.x) && (v1.y < v2.y) && (v1.z < v2.z);
}

inline bool operator <= (const FP3 & v1, const FP3 & v2)
{
    return (v1.x <= v2.x) && (v1.y <= v2.y) && (v1.z <= v2.z);
}

inline bool operator > (const FP3 & v1, const FP3 & v2)
{
    return (v1.x > v2.x) && (v1.y > v2.y) && (v1.z > v2.z);
}

inline bool operator >= (const FP3 & v1, const FP3 & v2)
{
    return (v1.x >= v2.x) && (v1.y >= v2.y) && (v1.z >= v2.z);
}

inline const Int3 truncate(const FP3& v)
{
    return Int3((int)v.x, (int)v.y, (int)v.z);
}

inline FP dot(const FP3& v1, const FP3& v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline FP SP(const FP3& v1, const FP3& v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

inline const FP3 cross(const FP3& v1, const FP3& v2)
{
    return FP3(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z,
        v1.x * v2.y - v1.y * v2.x);
}

inline const FP3 VP(const FP3& v1, const FP3& v2)
{
    return FP3(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z,
        v1.x * v2.y - v1.y * v2.x);
}

inline FP dist(const FP3& v1, const FP3& v2)
{
    return (v1 - v2).norm();
}

inline FP sqr(const FP3& v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

inline std::ostream& operator<<(std::ostream& out, const FP3& v)
{
    return out << "(" << v.x << ", " << v.y << ", " << v.z << ")";
}

} // namespace pica


#endif
