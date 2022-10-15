#pragma once

#include <algorithm>   

template<typename Type>
struct Vec2 {
    Vec2() = default;

    explicit Vec2(Type x) : x(x), y(x) {}

    Vec2(Type x, Type y) : x(x), y(y) {}

    Vec2 &operator+=(const Vec2 &b) {
        x += b.x;
        y += b.y;
        return *this;
    }

    Vec2 operator+(const Vec2 &b) const {
        Vec2 ret = *this;
        ret += b;
        return ret;
    }

    Vec2 &operator-=(const Vec2 &b) {
        x -= b.x;
        y -= b.x;
        return *this;
    }

    Vec2 operator-(const Vec2 &b) const {
        Vec2 ret = *this;
        ret -= b;
        return ret;
    }

    Vec2 operator-() const {
        return Vec2(-x, -y);
    }

    Vec2 operator*(double b) const {
        return Vec2(x * b, y * b);
    }

    Vec2 operator/(double b) const {
        return Vec2(x / b, y / b);
    }

    Vec2 operator/(Vec2 b) const {
        return Vec2(x / b.x, y / b.y);
    }

    Vec2 operator*(const Vec2 &b) const {
        return Vec2(x * b.x, y * b.y);
    }

    Type x, y;
};

template<typename Type>
struct Vec3 {
    Vec3() = default;

    explicit Vec3(Type x) : x(x), y(x), z(x) {}

    Vec3(Type x, Type y, Type z) : x(x), y(y), z(z) {}

    Vec3 operator+(const Vec3 &b) const {
        return Vec3(x + b.x, y + b.y, z + b.z);
    }

    Vec3 operator-(const Vec3 &b) const {
        return Vec3(x - b.x, y - b.y, z - b.z);
    }

    Vec3 operator-() const {
        return Vec3(-x, -y, -z);
    }

    Vec3 operator*(double b) const {
        return Vec3(x * b, y * b, z * b);
    }

    Vec3 operator/(double b) const {
        return Vec3(x / b, y / b, z / b);
    }

    Vec3 operator*(const Vec3 &b) const {
        return Vec3(x * b.x, y * b.y, z * b.z);
    }

    Vec3 &operator+=(const Vec3 &b) {
        x += b.x;
        y += b.y;
        z += b.z;
        return *this;
    }

    Type x, y, z;
};

template<typename Type>
struct Vec4 {
    Vec4() = default;

    explicit Vec4(Type x) : x(x), y(x), z(x), w(x) {}

    Vec4(Type x, Type y, Type z, Type w) : x(x), y(y), z(z), w(w) {}

    Type x, y, z, w;
};

template<typename Type>
Type dot(Vec3<Type> a, Vec3<Type> b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

template<typename Type>
Type length(Vec2<Type> a) {
    return std::sqrtf(a.x * a.x + a.y * a.y);
}

template<typename Type>
Type length(Vec3<Type> a) {
    return std::sqrtf(a.x * a.x + a.y * a.y + a.z * a.z);
}

template<typename Type>
Vec3<Type> cross(const Vec3<Type> &a, const Vec3<Type> &b) {
    return Vec3<Type>(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

template<typename Type>
Vec3<Type> normalize(Vec3<Type> a) {
    return a / (std::sqrt(a.x * a.x + a.y * a.y + a.z * a.z));
}

template<typename Type>
Vec3<Type> SphericalToVector(Type mu, Type phi) {
    return Vec3<Type>(std::sin(mu) * std::cos(phi),
                       std::sin(mu) * std::sin(phi),
                       std::cos(mu));
}

template<typename Type>
Vec3<Type> max(const Vec3<Type> &a, const Vec3<Type> &b) {
    return Vec3<Type>(std::max(a.x, b.x),
                       std::max(a.y, b.y),
                       std::max(a.z, b.z));
}

template<typename Type>
Vec3<Type> min(const Vec3<Type> &a, const Vec3<Type> &b) {
    return Vec3<Type>(std::min(a.x, b.x),
                       std::min(a.y, b.y),
                       std::min(a.z, b.z));
}

template<typename Type>
Vec3<Type> exp(const Vec3<Type> &vec) {
    return Vec3<Type>{ std::exp(vec.x), std::exp(vec.y) , std::exp(vec.z) };
}

using Vec2f = Vec2<float>;
using Vec3f = Vec3<float>;
using Vec4f = Vec4<float>;

using Vec2d = Vec2<double>;
using Vec3d = Vec3<double>;
using Vec4d = Vec4<double>;

