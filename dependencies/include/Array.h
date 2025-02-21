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
#include "SignalMath.h"

namespace Com
{

    template <typename NUMERIC> // INT, //DOUBLE //SHORT ... std::complex
    class Array : public std::vector<NUMERIC>
    {
    public:
        Array(std::initializer_list<NUMERIC> init);
        Array(const std::vector<NUMERIC> &vec);
        Array(std::size_t size);
        Array() = default;

        static Array<NUMERIC> arange(NUMERIC start, NUMERIC end, NUMERIC step);
        static Array<NUMERIC> linespace(NUMERIC start, NUMERIC end, size_t num);
        
        static Array<NUMERIC> ones(size_t N);

        virtual Array operator+(const Array &other) const;
        virtual Array operator-(const Array &other) const;
        virtual Array operator*(const Array &other) const;
        virtual Array operator/(const Array &other) const;

        virtual Array operator+(NUMERIC other) const;
        virtual Array operator-(NUMERIC other) const;
        virtual Array operator*(NUMERIC other) const;
        virtual Array operator/(NUMERIC other) const;

        Array slice(size_t start, size_t end, size_t step) const;
        Array operator()(size_t start, size_t end, size_t step) const;

        Array convolution(const Array &other) const;
        Array crossCorrelation(const Array &other) const;
        Array autoCorrelation() const;
        Array fourierTransform() const;
        Array applyFilter(const Array &filter) const;
        Array upSample(size_t factor) const;
        Array downSample(size_t factor) const;

        Array centeredDomain() const;

        double energy() const;
        double power() const;


        template <typename TARGET_TYPE>
        Array<TARGET_TYPE> convert() const;

        Array<NUMERIC> apply(std::function<NUMERIC(NUMERIC)> func) const;
        Array<NUMERIC> apply(std::function<NUMERIC(const NUMERIC &, double)> func, double x) const;
        Array<NUMERIC> apply(std::function<std::vector<NUMERIC>(NUMERIC, NUMERIC)> func, const std::vector<NUMERIC> &other) const;
        NUMERIC apply(std::function<NUMERIC(const std::vector<NUMERIC> &)> func) const;
        Array<double> apply(std::function<std::vector<double>(const std::vector<std::complex<double>> &)> func) const;
        Array<NUMERIC> apply(std::function<NUMERIC(NUMERIC, std::unordered_map<std::string, NUMERIC> &)> func,
                             std::unordered_map<std::string, NUMERIC> &params) const;

        std::vector<NUMERIC> toStdVector();
        std::string toString();

    private:
        class View
        {
        private:
            const Array &array;
            size_t _start;
            size_t _end;
            size_t step;

        public:
            View(const Array &Array, size_t Start, size_t End, size_t Step)
                : array(Array), _start(Start), _end(End), step(Step) {}

            class Iterator
            {
            public:
                Iterator(const Array &Array, size_t Index, size_t Step)
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
                size_t index;
                size_t step;
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
            std::vector<size_t> indices;

        public:
            SliceProxy(Array<NUMERIC> &parent, size_t start, size_t end, size_t step) : parent(parent)
            {
                if (step == 0)
                    throw std::invalid_argument("Step cannot be zero.");
                if (start >= parent.size() || end >= parent.size() || start > end)
                    throw std::out_of_range("Invalid slice indices.");

                for (size_t i = start; i < end; i += step)
                {
                    indices.push_back(i);
                }
            }

            SliceProxy &operator=(NUMERIC value)
            {
                for (size_t idx : indices)
                {
                    parent[idx] = value;
                }
                return *this;
            }
        };
        SliceProxy operator()(size_t start, size_t end, size_t step)
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

        for (size_t i = 0; i < num; ++i)
        {
            result [i] = (start + i * step);
        }
        return result;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::operator+(const Array &other) const
    {
        if (this->size() != other.size())
            throw std::invalid_argument("Vectors must have the same size");

        Array result(this->size());
        for (size_t i = 0; i < this->size(); ++i)
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
        for (size_t i = 0; i < this->size(); ++i)
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
        for (size_t i = 0; i < this->size(); ++i)
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
        for (size_t i = 0; i < this->size(); ++i)
        {
            if (std::abs(other[i]) == 0)
                throw std::invalid_argument("Division by zero");
            result[i] = this->at(i) / other[i];
        }
        return result;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::operator+(NUMERIC other) const
    {
        Array result(this->size());
        for (size_t i = 0; i < this->size(); ++i)
        {
            result[i] = this->at(i) + other;
        }
        return result;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::operator-(NUMERIC other) const
    {
        Array result(this->size());
        for (size_t i = 0; i < this->size(); ++i)
        {
            result[i] = this->at(i) - other;
        }
        return result;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::operator*(NUMERIC other) const
    {
        Array result(this->size());
        for (size_t i = 0; i < this->size(); ++i)
        {
            result[i] = this->at(i) * other;
        }
        return result;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::operator/(NUMERIC other) const
    {
        if (std::abs(other) == 0)
            throw std::invalid_argument("Division by zero");

        Array result(this->size());
        for (size_t i = 0; i < this->size(); ++i)
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
        size_t len1 = this->size();
        size_t len2 = other.size();

        std::vector<NUMERIC> result(len1 + len2 - 1, 0);

        for (size_t i = 0; i < len1; ++i)
        {
            for (size_t j = 0; j < len2; ++j)
            {
                result[i + j] += (*this)[i] * other[j];
            }
        }

        return result;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::crossCorrelation(const Array &other) const
    {
        size_t len1 = this->size();
        size_t len2 = other.size();

        std::vector<NUMERIC> result(len1 + len2 - 1, 0);

        for (size_t m = 0; m < result.size(); ++m)
        {
            for (size_t n = 0; n < len1; ++n)
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
    inline Array<NUMERIC> Array<NUMERIC>::fourierTransform() const
    {
        return Array();
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::applyFilter(const Array &filter) const
    {
        return this->convolution(filter);
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::upSample(size_t factor) const
    {
        if (factor < 1)
        {
            throw std::invalid_argument("Upsampling factor must be at least 1.");
        }

        std::vector<NUMERIC> result;
        result.reserve(this->size() * factor);

        for (size_t i = 0; i < this->size(); ++i)
        {
            result.push_back((*this)[i]); // Original sample
            for (size_t j = 1; j < factor; ++j)
            {
                result.push_back(0); // Insert zeros
            }
        }
        return result;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::downSample(size_t factor) const
    {
        if (factor < 1)
        {
            throw std::invalid_argument("Downsampling factor must be at least 1.");
        }

        std::vector<NUMERIC> result;
        result.reserve(this->size() / factor);

        for (size_t i = 0; i < this->size(); i += factor)
        {
            result.push_back((*this)[i]); // Keep every `factor`-th sample
        }

        return result;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::centeredDomain() const
    {
        double n = this->size();
        Array<double> domain = Array<double>::linespace(-n/2 , n/2, n);
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
    inline Array<NUMERIC> Array<NUMERIC>::ones(size_t N)
    {
        Array<NUMERIC> vec (N);
        vec.resize(N);
        for (size_t i = 0; i < N; i++)
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
        for (size_t i = 0; i < data.size() && i < other.size(); ++i)
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
    inline Array<double> Array<NUMERIC>::apply(std::function<std::vector<double>(const std::vector<std::complex<double>> &)> func) const
    {
        std::vector<double> result = func(*this);
        return Array<double>(result);
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::apply(std::function<NUMERIC(NUMERIC, std::unordered_map<std::string, NUMERIC> &)> func,
                                                std::unordered_map<std::string, NUMERIC> &params) const
    {
        std::vector<NUMERIC> result;
        for (const auto &elem : *this)
        {
            result.push_back(func(static_cast<NUMERIC>(elem), params));
        }
        return Array<NUMERIC>(result);
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
        for (size_t i = 0; i < this->size(); i++)
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
