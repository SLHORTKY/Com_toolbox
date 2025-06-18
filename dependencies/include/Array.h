#pragma once
#include <iostream>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <type_traits>
#include <functional>
#include <sstream>
#include <algorithm>
#include <complex>
#include <initializer_list>
#include <random>
#include <type_traits>
#include "SignalMath.h"

namespace Com
{
    template <typename T>
    struct is_complex : std::false_type
    {
    };

    template <typename T>
    struct is_complex<std::complex<T>> : std::true_type
    {
    };

    template <typename NUMERIC> // INT, //DOUBLE //SHORT ... std::complex
    class Array : public std::vector<NUMERIC>
    {
    public:
        Array(std::initializer_list<NUMERIC> init);
        Array(const std::vector<NUMERIC> &vec);
        Array(std::size_t size);
        Array() = default;

        static Array<NUMERIC> arange(NUMERIC start, NUMERIC end, NUMERIC step);
        static Array<NUMERIC> linespace(NUMERIC start, NUMERIC end, std::size_t num);

        static Array<NUMERIC> ones(std::size_t N);

        virtual Array operator+(const Array &other) const;
        virtual Array operator-(const Array &other) const;
        virtual Array operator*(const Array &other) const;
        virtual Array operator/(const Array &other) const;

        virtual Array operator+(NUMERIC other) const;
        virtual Array operator-(NUMERIC other) const;
        virtual Array operator*(NUMERIC other) const;
        virtual Array operator/(NUMERIC other) const;

        Array slice(std::size_t start, std::size_t end, std::size_t step) const;
        Array operator()(std::size_t start, std::size_t end, std::size_t step) const;

        Array convolution(const Array &other) const;
        Array crossCorrelation(const Array &other) const;
        Array autoCorrelation() const;
        Array applyFilter(const Array &filter) const;
        Array upSample(std::size_t factor) const;
        Array downSample(std::size_t factor) const;

        Array centeredDomain() const;

        double energy() const;
        double power() const;

        Array<double> real() const;
        Array<double> imag() const;

        Array<std::complex<double>> fft() const;
        Array<std::complex<double>> ifft() const;
        Array<std::complex<double>> fftshift() const;

        template <typename TARGET_TYPE>
        Array<TARGET_TYPE> convert() const;

        Array<NUMERIC> apply(std::function<NUMERIC(NUMERIC)> func) const;
        Array<NUMERIC> apply(std::function<NUMERIC(const NUMERIC &, double)> func, double x) const;
        Array<NUMERIC> apply(std::function<std::vector<NUMERIC>(NUMERIC, NUMERIC)> func, const std::vector<NUMERIC> &other) const;
        NUMERIC apply(std::function<NUMERIC(const std::vector<NUMERIC> &)> func) const;

        std::vector<NUMERIC> toStdVector();
        std::string toString();

    private:
        class View
        {
        private:
            const Array &array;
            std::size_t _start;
            std::size_t _end;
            std::size_t step;

        public:
            View(const Array &Array, std::size_t Start, std::size_t End, std::size_t Step)
                : array(Array), _start(Start), _end(End), step(Step) {}

            class Iterator
            {
            public:
                Iterator(const Array &Array, std::size_t Index, std::size_t Step)
                    : array(Array), index(Index), step(Step) {}

                Iterator &operator++()
                {
                    index += step;
                    return *this;
                }

                NUMERIC operator*() const
                {
                    return array[index];
                }

                bool operator!=(const Iterator &other) const
                {
                    return index != other.index;
                }

            private:
                const Array &array;
                std::size_t index;
                std::size_t step;
            };

            Iterator begin() const
            {
                return Iterator(array, _start, step);
            }

            Iterator end() const
            {
                return Iterator(array, _end, step);
            }
        };

    public:
        class SliceProxy
        {
        private:
            Array<NUMERIC> &parent;
            std::vector<std::size_t> indices;

        public:
            SliceProxy(Array<NUMERIC> &parent, std::size_t start, std::size_t end, std::size_t step) : parent(parent)
            {
                if (step == 0)
                    throw std::invalid_argument("Step cannot be zero.");
                if (start >= parent.size() || end >= parent.size() || start > end)
                    throw std::out_of_range("Invalid slice indices.");

                for (std::size_t i = start; i < end; i += step)
                {
                    indices.push_back(i);
                }
            }

            SliceProxy &operator=(NUMERIC value)
            {
                for (std::size_t idx : indices)
                {
                    parent[idx] = value;
                }
                return *this;
            }
        };
        SliceProxy operator()(std::size_t start, std::size_t end, std::size_t step)
        {
            return SliceProxy(*this, start, end, step);
        }
    };

    template <typename NUMERIC>
    inline Array<NUMERIC>::Array(std::initializer_list<NUMERIC> init)
        : std::vector<NUMERIC>(init) {}

    template <typename NUMERIC>
    inline Array<NUMERIC>::Array(const std::vector<NUMERIC> &vec)
        : std::vector<NUMERIC>(vec) {}

    template <typename NUMERIC>
    inline Array<NUMERIC>::Array(std::size_t size)
        : std::vector<NUMERIC>(size) {}

    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::arange(NUMERIC start, NUMERIC end, NUMERIC step)
    {
        if (step == 0)
            throw std::invalid_argument("Step cannot be zero");

        Array result;
        for (NUMERIC value = start; (step > 0 ? value < end : value > end); value += step)
        {
            result.push_back(value);
        }
        return result;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::linespace(NUMERIC start, NUMERIC end, std::size_t num)
    {
        if (num == 0)
            throw std::invalid_argument("Number of elements cannot be zero");

        Array<double> result = Array<double>::ones(num);
        NUMERIC step = (end - start) / static_cast<double>(num - 1);

        for (std::size_t i = 0; i < num; ++i)
        {
            result[i] = (start + i * step);
        }
        return result;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::operator+(const Array &other) const
    {
        if (this->size() != other.size())
            throw std::invalid_argument("Vectors must have the same size");

        Array result(this->size());
        for (std::size_t i = 0; i < this->size(); ++i)
        {
            result[i] = this->at(i) + other[i];
        }
        return result;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::operator-(const Array &other) const
    {
        if (this->size() != other.size())
            throw std::invalid_argument("Vectors must have the same size");

        Array result(this->size());
        for (std::size_t i = 0; i < this->size(); ++i)
        {
            result[i] = this->at(i) - other[i];
        }
        return result;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::operator*(const Array &other) const
    {
        if (this->size() != other.size())
            throw std::invalid_argument("Vectors must have the same size");

        Array result(this->size());
        for (std::size_t i = 0; i < this->size(); ++i)
        {
            result[i] = this->at(i) * other[i];
        }
        return result;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::operator/(const Array &other) const
    {
        if (this->size() != other.size())
            throw std::invalid_argument("Vectors must have the same size");

        Array result(this->size());
        for (std::size_t i = 0; i < this->size(); ++i)
        {

            if constexpr (is_complex<NUMERIC>::value)
            {
                if (std::abs(other[i]) == 0)
                    throw std::invalid_argument("Division by zero");
            }
            else
            {
                if (other[i] == 0)
                    throw std::invalid_argument("Division by zero");
            }

            result[i] = this->at(i) / other[i];
        }
        return result;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::operator+(NUMERIC other) const
    {
        Array result(this->size());
        for (std::size_t i = 0; i < this->size(); ++i)
        {
            result[i] = this->at(i) + other;
        }
        return result;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::operator-(NUMERIC other) const
    {
        Array result(this->size());
        for (std::size_t i = 0; i < this->size(); ++i)
        {
            result[i] = this->at(i) - other;
        }
        return result;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::operator*(NUMERIC other) const
    {
        Array result(this->size());
        for (std::size_t i = 0; i < this->size(); ++i)
        {
            result[i] = this->at(i) * other;
        }
        return result;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::operator/(NUMERIC other) const
    {

        if constexpr (is_complex<NUMERIC>::value)
        {
            if (std::abs(other) == 0)
                throw std::invalid_argument("Division by zero");
        }
        else
        {
            if (other == 0)
                throw std::invalid_argument("Division by zero");
        }
        Array result(this->size());
        for (std::size_t i = 0; i < this->size(); ++i)
        {
            result[i] = this->at(i) / other;
        }
        return result;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::operator()(std::size_t start, std::size_t end, std::size_t step) const
    {
        return this->slice(start, end, step);
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::convolution(const Array &other) const
    {
        std::size_t len1 = this->size();
        std::size_t len2 = other.size();

        std::vector<NUMERIC> result(len1 + len2 - 1, 0);

        for (std::size_t i = 0; i < len1; ++i)
        {
            for (std::size_t j = 0; j < len2; ++j)
            {
                result[i + j] += (*this)[i] * other[j];
            }
        }

        return result;
    }

    template <typename T>
    Array<T> operator+(T lhs, const Array<T> &rhs);

    template <typename T>
    Array<T> operator-(T lhs, const Array<T> &rhs);

    template <typename T>
    Array<T> operator*(T lhs, const Array<T> &rhs);

    template <typename T>
    Array<T> operator/(T lhs, const Array<T> &rhs);

    template <typename T>
    Array<T> operator+(T lhs, const Array<T> &rhs)
    {
        return rhs + lhs; // reuse the member operator+
    }

    template <typename T>
    Array<T> operator-(T lhs, const Array<T> &rhs)
    {
        // Assuming you want lhs - each element of rhs
        Array result(rhs.size());
        for (std::size_t i = 0; i < rhs.size(); ++i)
        {
            result[i] = lhs - rhs[i];
        }
        return result;
    }

    template <typename T>
    Array<T> operator*(T lhs, const Array<T> &rhs)
    {
        return rhs * lhs; // reuse the member operator*
    }

    template <typename T>
    Array<T> operator/(T lhs, const Array<T> &rhs)
    {
        Array result(rhs.size());
        for (std::size_t i = 0; i < rhs.size(); ++i)
        {
            result[i] = lhs / rhs[i];
        }
        return result;
    }

    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::crossCorrelation(const Array &other) const
    {
        std::size_t len1 = this->size();
        std::size_t len2 = other.size();

        std::vector<NUMERIC> result(len1 + len2 - 1, 0);

        for (std::size_t m = 0; m < result.size(); ++m)
        {
            for (std::size_t n = 0; n < len1; ++n)
            {
                if ((m - n) >= 0 && (m - n) < len2)
                {
                    result[m] += (*this)[n] * other[m - n];
                }
            }
        }

        return result;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::autoCorrelation() const
    {
        return this->crossCorrelation(*this);
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::applyFilter(const Array &filter) const
    {
        return this->convolution(filter);
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::upSample(std::size_t factor) const
    {
        if (factor < 1)
        {
            throw std::invalid_argument("Upsampling factor must be at least 1.");
        }

        std::vector<NUMERIC> result;
        result.reserve(this->size() * factor);

        for (std::size_t i = 0; i < this->size(); ++i)
        {
            result.push_back((*this)[i]); // Original sample
            for (std::size_t j = 1; j < factor; ++j)
            {
                result.push_back(0); // Insert zeros
            }
        }
        return result;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::downSample(std::size_t factor) const
    {
        if (factor < 1)
        {
            throw std::invalid_argument("Downsampling factor must be at least 1.");
        }

        std::vector<NUMERIC> result;
        result.reserve(this->size() / factor);

        for (std::size_t i = 0; i < this->size(); i += factor)
        {
            result.push_back((*this)[i]); // Keep every `factor`-th sample
        }

        return result;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::centeredDomain() const
    {
        double n = this->size();
        Array<double> domain = Array<double>::linespace(-n / 2, n / 2, n);
        return domain;
    }
    template <typename NUMERIC>
    inline double Array<NUMERIC>::energy() const
    {
        return this->apply(SignalMath::pow, 2.0).apply(SignalMath::sum);
    }
    template <typename NUMERIC>
    inline double Array<NUMERIC>::power() const
    {
        return this->apply(SignalMath::pow, 2.0).apply(SignalMath::mean);
    }
    template <typename NUMERIC>
    inline Array<std::complex<double>> Array<NUMERIC>::fft() const
    {
        if constexpr (std::is_same_v<NUMERIC, std::complex<double>>)
        {
            return SignalMath::fft_comp(*this);
        }
        else
        {
            return SignalMath::fft(*this);
        }
    }

    template <typename NUMERIC>
    inline Array<std::complex<double>> Array<NUMERIC>::ifft() const
    {
        return SignalMath::ifft(*this);
    }
    template <typename NUMERIC>
    inline Array<std::complex<double>> Array<NUMERIC>::fftshift() const
    {
        // Get the size of the array
        size_t N = this->size();

        // Output array with the same size
        Array<std::complex<double>> shifted(N);

        // Compute the midpoint (half the length)
        size_t mid = N / 2;

        // Shift the second half to the front
        for (size_t i = 0; i < N - mid; ++i)
            shifted[i] = std::complex<double>((*this)[mid + i]);

        // Shift the first half to the back
        for (size_t i = 0; i < mid; ++i)
            shifted[N - mid + i] = std::complex<double>((*this)[i]);

        return shifted;
    }
    template <typename NUMERIC>
    inline Array<double> Array<NUMERIC>::real() const
    {
        return SignalMath::real(*this);
    }
    template <typename NUMERIC>
    inline Array<double> Array<NUMERIC>::imag() const
    {
        return SignalMath::imag(*this);
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::ones(std::size_t N)
    {
        Array<NUMERIC> vec(N);
        vec.resize(N);
        for (std::size_t i = 0; i < N; i++)
        {
            vec[i] = static_cast<NUMERIC>(1);
        }
        return vec;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::slice(std::size_t start, std::size_t end, std::size_t step) const
    {
        {
            if (step == 0)
                throw std::invalid_argument("Step size cannot be zero");
            if (start >= this->size() || end > this->size() || start > end)
                throw std::out_of_range("Invalid slice indices");

            View view(*this, start, end, step);
            Array result;
            for (auto value : view)
            {
                result.push_back(value);
            }
            return result;
        }
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::apply(std::function<NUMERIC(NUMERIC)> func) const
    {
        std::vector<NUMERIC> result;
        for (const auto &elem : *this)
        {
            result.push_back(func(static_cast<NUMERIC>(elem)));
        }
        return Array<NUMERIC>(result);
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::apply(std::function<NUMERIC(const NUMERIC &, double)> func, double x) const
    {
        std::vector<NUMERIC> result;
        for (const auto &elem : *this)
        {
            result.push_back(func(elem, x));
        }
        return Array<NUMERIC>(result);
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::apply(std::function<std::vector<NUMERIC>(NUMERIC, NUMERIC)> func, const std::vector<NUMERIC> &other) const
    {
        std::vector<NUMERIC> result;
        std::vector<NUMERIC> data = *this;
        for (std::size_t i = 0; i < data.size() && i < other.size(); ++i)
        {
            auto temp_result = func(data[i], other[i]);
            result.insert(result.end(), temp_result.begin(), temp_result.end());
        }
        return Array<NUMERIC>(result);
    }
    template <typename NUMERIC>
    inline NUMERIC Array<NUMERIC>::apply(std::function<NUMERIC(const std::vector<NUMERIC> &)> func) const
    {
        return func(*this);
    }

    template <typename NUMERIC>
    inline std::vector<NUMERIC> Array<NUMERIC>::toStdVector()
    {
        std::vector<NUMERIC> vec = *this;
        return vec;
    }
    template <typename NUMERIC>
    inline std::string Array<NUMERIC>::toString()
    {
        std::stringstream ss;
        ss << "[";
        for (std::size_t i = 0; i < this->size(); i++)
        {
            ss << this->at(i);
            i < this->size() - 1 ? ss << ", " : ss << "";
        }
        ss << "]";
        return ss.str();
    }
    template <typename NUMERIC>
    template <typename TARGET_TYPE>
    inline Array<TARGET_TYPE> Array<NUMERIC>::convert() const
    {
        Array<TARGET_TYPE> result;
        result.reserve(this->size());

        for (const auto &value : *this)
        {
            result.push_back(static_cast<TARGET_TYPE>(value));
        }

        return result;
    }
} // namespace Com
