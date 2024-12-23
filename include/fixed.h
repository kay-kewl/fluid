#pragma once

#include <cstdint>
#include <limits>
#include <ostream>
#include <istream>
#include <array>
#include <algorithm>

template <size_t N, size_t K>
struct Fixed {
    using IntType = int64_t;
    static constexpr size_t Bits = N;
    static constexpr size_t Fraction = K;
    
    constexpr Fixed() noexcept : v(0) {}
    constexpr explicit Fixed(IntType v_ = 0) : v(v_) {}
    constexpr explicit Fixed(float f) : v(static_cast<IntType>(f * (1ULL << K))) {}
    constexpr explicit Fixed(double f) : v(static_cast<IntType>(f * (1ULL << K))) {}

    static constexpr Fixed from_raw(IntType x) {
        Fixed ret;
        ret.v = x;
        return ret;
    }

    IntType v;

    auto operator<=>(const Fixed&) const = default;
    bool operator==(const Fixed&) const = default;

    friend std::ostream& operator<<(std::ostream& out, const Fixed& f) {
        out << (static_cast<double>(f.v) / (1ULL << K));
        return out;
    }
    
    friend std::istream& operator>>(std::istream& in, Fixed& f) {
        double temp;
        in >> temp;
        f = Fixed(temp);
        return in;
    }
};

template <size_t N, size_t K>
Fixed<N, K> operator+(Fixed<N, K> a, Fixed<N, K> b) {
    return Fixed<N, K>::from_raw(a.v + b.v);
}

template <size_t N, size_t K>
Fixed<N, K> operator-(Fixed<N, K> a, Fixed<N, K> b) {
    return Fixed<N, K>::from_raw(a.v - b.v);
}

template <size_t N, size_t K>
Fixed<N, K> operator*(Fixed<N, K> a, Fixed<N, K> b) {
    return Fixed<N, K>::from_raw((a.v * b.v) >> K);
}

template <size_t N, size_t K>
Fixed<N, K> operator/(Fixed<N, K> a, Fixed<N, K> b) {
    return Fixed<N, K>::from_raw((a.v << K) / b.v);
}

template <size_t N, size_t K>
Fixed<N, K>& operator+=(Fixed<N, K>& a, Fixed<N, K> b) {
    a = a + b;
    return a;
}

template <size_t N, size_t K>
Fixed<N, K>& operator-=(Fixed<N, K>& a, Fixed<N, K> b) {
    a = a - b;
    return a;
}

template <size_t N, size_t K>
Fixed<N, K>& operator*=(Fixed<N, K>& a, Fixed<N, K> b) {
    a = a * b;
    return a;
}

template <size_t N, size_t K>
Fixed<N, K>& operator/=(Fixed<N, K>& a, Fixed<N, K> b) {
    a = a / b;
    return a;
}

template <size_t N, size_t K>
Fixed<N, K> operator-(Fixed<N, K> x) {
    return Fixed<N, K>::from_raw(-x.v);
}

template <size_t N, size_t K>
Fixed<N, K> abs(Fixed<N, K> x) {
    if (x.v < 0) {
        x.v = -x.v;
    }
    return x;
}

template <size_t N, size_t K>
Fixed<N, K> operator*(Fixed<N, K> a, double b) {
    return Fixed<N, K>(static_cast<double>(a.v) * b / (1ULL << K));
}

template <size_t N, size_t K>
Fixed<N, K> operator*(double a, Fixed<N, K> b) {
    return Fixed<N, K>(a * static_cast<double>(b.v) / (1ULL << K));
}

template <size_t N, size_t K>
Fixed<N, K> operator/(Fixed<N, K> a, double b) {
    return Fixed<N, K>(static_cast<double>(a.v) / ((1ULL << K) * b));
}

template <size_t N, size_t K>
Fixed<N, K>& operator*=(Fixed<N, K>& a, double b) {
    a = a * b;
    return a;
}

template <size_t N, size_t K>
Fixed<N, K>& operator/=(Fixed<N, K>& a, double b) {
    a = a / b;
    return a;
}

template <size_t N, size_t K>
Fixed<N, K> operator*(Fixed<N, K> a, float b) {
    return a * static_cast<double>(b);
}

template <size_t N, size_t K>
Fixed<N, K> operator*(float a, Fixed<N, K> b) {
    return static_cast<double>(a) * b;
}

template <size_t N, size_t K>
Fixed<N, K> operator/(Fixed<N, K> a, float b) {
    return a / static_cast<double>(b);
}

template <size_t N, size_t K>
Fixed<N, K>& operator*=(Fixed<N, K>& a, float b) {
    a = a * b;
    return a;
}

template <size_t N, size_t K>
Fixed<N, K>& operator/=(Fixed<N, K>& a, float b) {
    a = a / b;
    return a;
}

template<size_t N, size_t K>
struct FastFixed : public Fixed<N,K> {
    using Base = Fixed<N,K>;
    using typename Base::IntType;
    
    static constexpr size_t Bits = N;
    static constexpr size_t Fraction = K;
    
    using Base::Base;
    using Base::operator=;
    using Base::from_raw;
    
    friend std::ostream& operator<<(std::ostream& out, const FastFixed& f) {
        return out << static_cast<const Base&>(f);
    }
    
    friend std::istream& operator>>(std::istream& in, FastFixed& f) {
        Base& base = f;
        return in >> base;
    }
};