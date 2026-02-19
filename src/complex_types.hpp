#pragma once

// These are included to group all compatible FFT input types in this header
#include <complex>
#include <arm_neon.h>

namespace FFT {
    // Prefered Complex Type for fastest results, as it has best alignment without intrinsics.
    //  This gives the compiler the greatest freedom in selecting SIMD instructions.
    struct alignas(8) Complex {
        float data[2];
        
        // Accessors
        float& operator[](size_t i) { return data[i]; }
        const float& operator[](size_t i) const { return data[i]; }
        
        // Default constructor
        constexpr Complex() : data{0, 0} {}
        
        // Initialization constructor
        constexpr Complex(float r, float i = 0.0f) : data{r, i} {}
        
        // Operators
        friend constexpr Complex operator+(const Complex& a, const Complex& b) noexcept {
            return {a[0] + b[0], a[1] + b[1]};
        }
        
        friend constexpr Complex operator-(const Complex& a, const Complex& b) noexcept {
            return {a[0] - b[0], a[1] - b[1]};
        }

        friend constexpr Complex operator*(const Complex& a, float b) noexcept {
            return {a[0] * b, a[1] * b};
        }
        
        friend constexpr Complex operator/(const Complex& a, float b) noexcept {
            return {a[0] / b, a[1] / b};
        }
        
        friend constexpr Complex operator*(const Complex& a, const Complex& b) noexcept {
            return {a[0] * b[0], a[1] * b[1]};
        }

        Complex& operator+=(const Complex& b) noexcept {
            data[0] += b[0];
            data[1] += b[1];
            return *this;
        }
        
        Complex& operator/=(float b) noexcept {
            data[0] /= b;
            data[1] /= b;
            return *this;
        }
        
        // Implicit conversion from std::complex for compatibility
        Complex(const std::complex<float>& c) : data{c.real(), c.imag()} {}
        operator std::complex<float>() const { return {data[0], data[1]}; }
    };

    // This is a wrapper type for std::complex to allow it to be used as an input to FFTPlan.
    template <typename T>
    struct StdComplexWrap : std::complex<T> {
        // Conversion to std::complex
        operator std::complex<T>&() { return *this; }
        operator const std::complex<T>&() const { return *this; }

        // Construct from std::complex
        constexpr StdComplexWrap(const std::complex<T>& value = {}) : std::complex<T>(value) {}
        constexpr StdComplexWrap(T r, T i = 0.0f) : std::complex<T>(r, i) {}
        
        // Accessor for FFTPlan compatability
        T& operator[](size_t i) { return reinterpret_cast<T*>(this)[i]; }
        const T& operator[](size_t i) const { return reinterpret_cast<const T*>(this)[i]; }
    };
}