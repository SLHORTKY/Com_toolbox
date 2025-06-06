#pragma once
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cmath>

namespace Com
{
    class SignalMath
    {
    public:
        static double sum(const std::vector<double> &vec);
        static double mean(const std::vector<double> &vec);
        static double median(const std::vector<double> &vec);
        static double max(const std::vector<double> &vec);
        static double min(const std::vector<double> &vec);

        static std::vector<double> real(const std::vector<std::complex<double>> &vec);
        static std::vector<double> imag(const std::vector<std::complex<double>> &vec);

        static std::complex<double> rotate(const std::complex<double> &x, double angle);

        static std::vector<std::complex<double>> fft(const std::vector<std::complex<double>> &input);
        static std::vector<std::complex<double>> ifft(const std::vector<std::complex<double>> &input);

        static double deg2rad(double x);
        static double rad2deg(double x);
        static double normalize(double x, double value); // normalize relative to a value;

        static double sqrt(double x);               // Square root of x
        static double abs(double x);                // Absolute value of x
        static double ceil(double x);               // Ceiling function (round up)
        static double floor(double x);              // Floor function (round down)
        static double fmod(double x, double y);     // Remainder of division (x % y)
        static double pow(double base, double exp); // x raised to the power of y
        static double exp(double x);                // e raised to the power of x (exp(x))
        static double log(double x);                // Natural logarithm (ln(x))
        static double log10(double x);              // Base-10 logarithm (log(x))
        static double sin(double x);                // Sine of x (in radians)
        static double cos(double x);                // Cosine of x (in radians)
        static double tan(double x);                // Tangent of x (in radians)
        static double asin(double x);               // Arc sine (inverse of sin)
        static double acos(double x);               // Arc cosine (inverse of cos)
        static double atan(double x);               // Arc tangent (inverse of tan)
        static double atan2(double y, double x);    // Arc tangent (imag, real)
        static double sinh(double x);               // Hyperbolic Sin
        static double cosh(double x);               // Hyperbolic Cos
        static double tanh(double x);               // Hyperbolic Tan
        static double hypot(double x, double y);    // Hypotenuse
    };

    inline double SignalMath::sum(const std::vector<double> &vec)
    {
        double result = 0.0;

        for (size_t i = 0; i < vec.size(); i++)
        {
            result += vec.at(i);
        }
        return result;
    }

    inline double SignalMath::mean(const std::vector<double> &vec)
    {
        return sum(vec) / vec.size();
    }

    inline double SignalMath::median(const std::vector<double> &vec)
    {
        if (vec.empty())
        {
            throw std::invalid_argument("Vector is empty");
        }

        std::vector<double> sorted_vec = vec; // Copy the vector to preserve the original order
        std::sort(sorted_vec.begin(), sorted_vec.end());

        size_t size = sorted_vec.size();
        if (size % 2 == 0)
        {
            // If even number of elements, return the average of the two middle values
            return (sorted_vec[size / 2 - 1] + sorted_vec[size / 2]) / 2.0;
        }
        else
        {
            // If odd number of elements, return the middle value
            return sorted_vec[size / 2];
        }
    }

    inline double SignalMath::max(const std::vector<double> &vec)
    {
        if (vec.empty())
        {
            throw std::invalid_argument("Vector is empty");
        }

        return *std::max_element(vec.begin(), vec.end());
    }

    inline double SignalMath::min(const std::vector<double> &vec)
    {
        if (vec.empty())
        {
            throw std::invalid_argument("Vector is empty");
        }

        return *std::min_element(vec.begin(), vec.end());
    }

    inline std::vector<double> SignalMath::real(const std::vector<std::complex<double>> &vec)
    {
        std::vector<double> real_vector;
        if (vec.empty())
        {
            throw std::invalid_argument("Vector is empty");
        }
        for (size_t i = 0; i < vec.size(); i++)
        {
            real_vector.push_back(vec[i].real());
        }
        return real_vector;
    }

    inline std::vector<double> SignalMath::imag(const std::vector<std::complex<double>> &vec)
    {
        std::vector<double> real_vector;
        if (vec.empty())
        {
            throw std::invalid_argument("Vector is empty");
        }
        for (size_t i = 0; i < vec.size(); i++)
        {
            real_vector.push_back(vec[i].imag());
        }
        return real_vector;
    }

    inline std::complex<double> SignalMath::rotate(const std::complex<double> &x, double angle)
    {
        double rad = SignalMath::deg2rad(angle);
        return x * std::complex<double>(cos(rad), sin(rad));
    }

    inline std::vector<std::complex<double>> SignalMath::fft(const std::vector<std::complex<double>> &input)
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

        auto fft_even = fft(even);
        auto fft_odd = fft(odd);

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

        // Take conjugate of input
        std::vector<std::complex<double>> conjugated_input(N);
        for (std::size_t i = 0; i < N; ++i)
        {
            conjugated_input[i] = conj(input[i]);
        }

        // Compute forward FFT of conjugated input
        auto fft_result = fft(conjugated_input);

        // Take conjugate again and scale by 1/N
        std::vector<std::complex<double>> result(N);
        for (std::size_t i = 0; i < N; ++i)
        {
            result[i] = std::conj(fft_result[i]) / static_cast<double>(N);
        }

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

    inline double SignalMath::sqrt(double x)
    {
        if (x < 0)
            throw std::invalid_argument("Cannot calculate the square root of a negative number.");
        return std::sqrt(x);
    }

    inline double SignalMath::abs(double x)
    {
        return std::abs(x);
    }

    inline double SignalMath::ceil(double x)
    {
        return std::ceil(x);
    }

    inline double SignalMath::floor(double x)
    {
        return std::floor(x);
    }

    inline double SignalMath::fmod(double x, double y)
    {
        if (y == 0)
            throw std::invalid_argument("Division by zero.");
        return std::fmod(x, y);
    }

    inline double SignalMath::pow(double base, double exp)
    {
        return std::pow(base, exp);
    }

    inline double SignalMath::exp(double x)
    {
        return std::exp(x);
    }

    inline double SignalMath::log(double x)
    {
        if (x <= 0)
            throw std::invalid_argument("Logarithm undefined for non-positive values.");
        return std::log(x);
    }

    inline double SignalMath::log10(double x)
    {
        if (x <= 0)
            throw std::invalid_argument("Logarithm undefined for non-positive values.");
        return std::log10(x);
    }

    inline double SignalMath::sin(double x)
    {
        return std::sin(x);
    }

    inline double SignalMath::cos(double x)
    {
        return std::cos(x);
    }

    inline double SignalMath::tan(double x)
    {
        return std::tan(x);
    }

    inline double SignalMath::asin(double x)
    {
        if (x < -1 || x > 1)
            throw std::invalid_argument("Input must be between -1 and 1.");
        return std::asin(x);
    }

    inline double SignalMath::acos(double x)
    {
        if (x < -1 || x > 1)
            throw std::invalid_argument("Input must be between -1 and 1.");
        return std::acos(x);
    }

    inline double SignalMath::atan(double x)
    {
        return std::atan(x);
    }

    inline double SignalMath::atan2(double y, double x)
    {
        return std::atan2(y, x);
    }

    inline double SignalMath::sinh(double x)
    {
        return std::sinh(x);
    }

    inline double SignalMath::cosh(double x)
    {
        return std::cosh(x);
    }

    inline double SignalMath::tanh(double x)
    {
        return std::tanh(x);
    }

    inline double SignalMath::hypot(double x, double y)
    {
        return std::hypot(x, y);
    }

} // namespace Com
