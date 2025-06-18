#pragma once
#include <vector>
#include <complex>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <type_traits>

namespace Com
{

    class SignalMath
    {
    public:
        template <typename T>
        static T sum(const std::vector<T> &vec);
        template <typename T>
        static T mean(const std::vector<T> &vec);
        template <typename T>
        static T max(const std::vector<T> &vec);
        template <typename T>
        static T min(const std::vector<T> &vec);

        template <typename T>
        static std::vector<double> real(const std::vector<T> &vec);
        template <typename T>
        static std::vector<double> imag(const std::vector<T> &vec);

        template <typename T>
        static std::vector<std::complex<double>> fft(const std::vector<T> &input);
        static std::vector<std::complex<double>> fft_comp(const std::vector<std::complex<double>> &input);
        static std::vector<std::complex<double>> ifft(const std::vector<std::complex<double>> &input);

        static std::complex<double> rotate(const std::complex<double> &x, double angle);

        static double deg2rad(double x);
        static double rad2deg(double x);
        static double normalize(double x, double value);

        template <typename T>
        static T sqrt(T x);
        template <typename T>
        static T abs(T x);
        template <typename T>
        static T ceil(T x);
        template <typename T>
        static T floor(T x);
        template <typename T>
        static T fmod(T x, T y);
        template <typename T>
        static T pow(T base, T exp);
        template <typename T>
        static T exp(T x);
        template <typename T>
        static T log(T x);
        template <typename T>
        static T log10(T x);

        template <typename T>
        static T sin(T x);
        template <typename T>
        static T cos(T x);
        template <typename T>
        static T tan(T x);
        template <typename T>
        static T asin(T x);
        template <typename T>
        static T acos(T x);
        template <typename T>
        static T atan(T x);
        static double atan2(double y, double x);

        template <typename T>
        static T sinh(T x);
        template <typename T>
        static T cosh(T x);
        template <typename T>
        static T tanh(T x);

        static double hypot(double x, double y);
    };


    template <>
    inline std::vector<double> SignalMath::real(const std::vector<std::complex<double>> &vec);
    template <>
    inline std::vector<double> SignalMath::imag(const std::vector<std::complex<double>> &vec);



    template <typename T>
    inline T SignalMath::sum(const std::vector<T> &vec)
    {
        T result = T(0);
        for (const auto &val : vec)
            result += val;
        return result;
    }

    template <typename T>
    inline T SignalMath::mean(const std::vector<T> &vec)
    {
        if (vec.empty())
            throw std::invalid_argument("Vector is empty");
        return sum(vec) / static_cast<double>(vec.size());
    }

    template <typename T>
    inline T SignalMath::max(const std::vector<T> &vec)
    {
        if (vec.empty())
            throw std::invalid_argument("Vector is empty");
        return *std::max_element(vec.begin(), vec.end());
    }

    template <typename T>
    inline T SignalMath::min(const std::vector<T> &vec)
    {
        if (vec.empty())
            throw std::invalid_argument("Vector is empty");
        return *std::min_element(vec.begin(), vec.end());
    }

    template <typename T>
    inline std::vector<double> SignalMath::real(const std::vector<T> &vec)
    {
        std::vector<double> result;
        result.reserve(vec.size());
        for (const auto &v : vec)
            result.push_back(static_cast<double>(v));
        return result;
    }

    template <>
    inline std::vector<double> SignalMath::real(const std::vector<std::complex<double>> &vec)
    {
        std::vector<double> result;
        result.reserve(vec.size());
        for (const auto &v : vec)
            result.push_back(v.real());
        return result;
    }

    template <typename T>
    inline std::vector<double> SignalMath::imag(const std::vector<T> &vec)
    {
        return std::vector<double>(vec.size(), 0.0);
    }

    template <>
    inline std::vector<double> SignalMath::imag(const std::vector<std::complex<double>> &vec)
    {
        std::vector<double> result;
        result.reserve(vec.size());
        for (const auto &v : vec)
            result.push_back(v.imag());
        return result;
    }

    inline std::complex<double> SignalMath::rotate(const std::complex<double> &x, double angle)
    {
        double rad = deg2rad(angle);
        return x * std::complex<double>(cos(rad), sin(rad));
    }

    template <typename T>
    inline std::vector<std::complex<double>> SignalMath::fft(const std::vector<T> &input)
    {
        std::vector<std::complex<double>> complex_input(input.begin(), input.end());
        return fft_comp(complex_input);
    }

    inline std::vector<std::complex<double>> SignalMath::fft_comp(const std::vector<std::complex<double>> &input)
    {
        std::size_t N = input.size();
        if (N <= 1)
            return input;

        std::vector<std::complex<double>> even(N / 2), odd(N / 2);
        for (std::size_t i = 0; i < N / 2; ++i)
        {
            even[i] = input[i * 2];
            odd[i] = input[i * 2 + 1];
        }

        auto fft_even = fft_comp(even);
        auto fft_odd = fft_comp(odd);

        std::vector<std::complex<double>> result(N);
        for (std::size_t k = 0; k < N / 2; ++k)
        {
            std::complex<double> t = std::polar(1.0, -2 * M_PI * k / N) * fft_odd[k];
            result[k] = fft_even[k] + t;
            result[k + N / 2] = fft_even[k] - t;
        }

        return result;
    }

    inline std::vector<std::complex<double>> SignalMath::ifft(const std::vector<std::complex<double>> &input)
    {
        std::size_t N = input.size();
        std::vector<std::complex<double>> conjugated_input(N);
        for (std::size_t i = 0; i < N; ++i)
            conjugated_input[i] = std::conj(input[i]);

        auto fft_result = fft_comp(conjugated_input);

        std::vector<std::complex<double>> result(N);
        for (std::size_t i = 0; i < N; ++i)
            result[i] = std::conj(fft_result[i]) / static_cast<double>(N);

        return result;
    }

    inline double SignalMath::deg2rad(double x)
    {
        return x / 180.0 * M_PI;
    }

    inline double SignalMath::rad2deg(double x)
    {
        return x * 180.0 / M_PI;
    }

    inline double SignalMath::normalize(double x, double value)
    {
        return x / value;
    }

    template <typename T>
    inline T SignalMath::sqrt(T x)
    {
        if constexpr (std::is_same_v<T, double>)
        {
            if (x < 0)
                throw std::invalid_argument("Cannot calculate the square root of a negative number.");
        }
        return std::sqrt(x);
    }

    template <typename T>
    inline T SignalMath::abs(T x)
    {
        return std::abs(x);
    }

    template <typename T>
    inline T SignalMath::ceil(T x)
    {
        if constexpr (std::is_same_v<T, std::complex<double>>)
            throw std::invalid_argument("ceil not defined for complex numbers");
        return std::ceil(x);
    }

    template <typename T>
    inline T SignalMath::floor(T x)
    {
        if constexpr (std::is_same_v<T, std::complex<double>>)
            throw std::invalid_argument("floor not defined for complex numbers");
        return std::floor(x);
    }

    template <typename T>
    inline T SignalMath::fmod(T x, T y)
    {
        if constexpr (std::is_same_v<T, std::complex<double>>)
            throw std::invalid_argument("fmod not defined for complex numbers");
        if (y == 0)
            throw std::invalid_argument("Division by zero.");
        return std::fmod(x, y);
    }

    template <typename T>
    inline T SignalMath::pow(T base, T exp)
    {
        return std::pow(base, exp);
    }

    template <typename T>
    inline T SignalMath::exp(T x)
    {
        return std::exp(x);
    }

    template <typename T>
    inline T SignalMath::log(T x)
    {
        if constexpr (std::is_same_v<T, double> && x <= 0)
            throw std::invalid_argument("Logarithm undefined for non-positive values.");
        return std::log(x);
    }

    template <typename T>
    inline T SignalMath::log10(T x)
    {
        if constexpr (std::is_same_v<T, double> && x <= 0)
            throw std::invalid_argument("Logarithm undefined for non-positive values.");
        return std::log10(x);
    }

    template <typename T>
    inline T SignalMath::sin(T x)
    {
        return std::sin(x);
    }

    template <typename T>
    inline T SignalMath::cos(T x)
    {
        return std::cos(x);
    }

    template <typename T>
    inline T SignalMath::tan(T x)
    {
        return std::tan(x);
    }

    template <typename T>
    inline T SignalMath::asin(T x)
    {
        if constexpr (std::is_same_v<T, double> && (x < -1 || x > 1))
            throw std::invalid_argument("Input must be between -1 and 1.");
        return std::asin(x);
    }

    template <typename T>
    inline T SignalMath::acos(T x)
    {
        if constexpr (std::is_same_v<T, double> && (x < -1 || x > 1))
            throw std::invalid_argument("Input must be between -1 and 1.");
        return std::acos(x);
    }

    template <typename T>
    inline T SignalMath::atan(T x)
    {
        return std::atan(x);
    }

    inline double SignalMath::atan2(double y, double x)
    {
        return std::atan2(y, x);
    }

    template <typename T>
    inline T SignalMath::sinh(T x)
    {
        return std::sinh(x);
    }

    template <typename T>
    inline T SignalMath::cosh(T x)
    {
        return std::cosh(x);
    }

    template <typename T>
    inline T SignalMath::tanh(T x)
    {
        return std::tanh(x);
    }

    inline double SignalMath::hypot(double x, double y)
    {
        return std::hypot(x, y);
    }

} // namespace Com
