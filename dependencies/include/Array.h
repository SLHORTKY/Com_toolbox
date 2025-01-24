#pragma once
#include <iostream>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <initializer_list>

namespace Com
{
    class Array : public std::vector<double>
    {

    public:
        Array(std::initializer_list<double> init);
        Array(const std::vector<double> &vec);
        Array(size_t size);
        Array() = default;

        static Array arange(double start, double end, double step);
        static Array linespace(double start, double end, size_t num);

        Array operator+(const Array &other) const;
        Array operator-(const Array &other) const;
        Array operator*(const Array &other) const;
        Array operator/(const Array &other) const;

        Array operator+(double other) const;
        Array operator-(double other) const;
        Array operator*(double other) const;
        Array operator/(double other) const;

        Array slice(size_t start, size_t end, size_t step) const;

        Array operator()(size_t start, size_t end, size_t step) const;

        Array apply(double (*func)(double)) const;
        Array apply(double (*func)(double, double), const Array &other) const;
        Array apply(double (*func)(double, std::unordered_map<std::string, double> &params), std::unordered_map<std::string, double> &params) const;

        std::vector<double> toStdVector();
        std::string toString();
    };

    inline Array::Array(std::initializer_list<double> init)
        : std::vector<double>(init) {}

    inline Array::Array(const std::vector<double> &vec)
        : std::vector<double>(vec) {}

    inline Array::Array(size_t size)
        : std::vector<double>(size) {}

    inline Array Array::arange(double start, double end, double step)
    {
        if (step == 0)
            throw std::invalid_argument("Step cannot be zero");

        Array result;
        for (double value = start; (step > 0 ? value < end : value > end); value += step)
        {
            result.push_back(value);
        }
        return result;
    }
    inline Array Array::linespace(double start, double end, size_t num)
    {
        if (num == 0)
            throw std::invalid_argument("Number of elements cannot be zero");

        Array result;
        double step = (end - start) / static_cast<double>(num - 1);
        for (size_t i = 0; i < num; ++i)
        {
            result.push_back(start + i * step);
        }
        return result;
    }
    inline Array Array::operator+(const Array &other) const
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
    inline Array Array::operator-(const Array &other) const
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
    inline Array Array::operator*(const Array &other) const
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
    inline Array Array::operator/(const Array &other) const
    {
        if (this->size() != other.size())
            throw std::invalid_argument("Vectors must have the same size");

        Array result(this->size());
        for (size_t i = 0; i < this->size(); ++i)
        {
            if (other[i] == 0)
                throw std::invalid_argument("Division by zero");
            result[i] = this->at(i) / other[i];
        }
        return result;
    }
    inline Array Array::operator+(double other) const
    {
        Array result(this->size());
        for (size_t i = 0; i < this->size(); ++i)
        {
            result[i] = this->at(i) + other;
        }
        return result;
    }
    inline Array Array::operator-(double other) const
    {
        Array result(this->size());
        for (size_t i = 0; i < this->size(); ++i)
        {
            result[i] = this->at(i) - other;
        }
        return result;
    }
    inline Array Array::operator*(double other) const
    {
        Array result(this->size());
        for (size_t i = 0; i < this->size(); ++i)
        {
            result[i] = this->at(i) * other;
        }
        return result;
    }
    inline Array Array::operator/(double other) const
    {
        if (other == 0)
            throw std::invalid_argument("Division by zero");

        Array result(this->size());
        for (size_t i = 0; i < this->size(); ++i)
        {
            result[i] = this->at(i) / other;
        }
        return result;
    }
    inline Array Array::operator()(size_t start, size_t end, size_t step) const
    {
        return this->slice(start, end, step);
    }
    inline Array Array::slice(size_t start, size_t end, size_t step = 1) const
    {
        if (start >= this->size() || end > this->size() || start > end || step == 0)
        {
            throw std::out_of_range("Invalid slice indices or step size");
        }

        Array result;
        for (size_t i = start; i < end; i += step)
        {
            result.push_back(this->at(i));
        }
        return result;
    }
    inline Array Array::apply(double (*func)(double)) const
    {
        Array vec;
        for (size_t i = 0; i < this->size(); i++)
        {
            vec.push_back(func(this->at(i)));
        }
        return vec;
    }
    inline Array Array::apply(double (*func)(double, double), const Array &other) const
    {
        Array vec;
        for (size_t i = 0; i < this->size(); i++)
        {
            vec.push_back(func(this->at(i), other.at(i)));
        }
        return vec;
    }
    inline Array Array::apply(double (*func)(double, std::unordered_map<std::string, double> &params), std::unordered_map<std::string, double> &params) const
    {
        std::vector<double> vec;
        for (size_t i = 0; i < this->size(); i++)
        {
            vec.push_back(func(this->at(i), params));
        }
        return vec;
    }
    inline std::vector<double> Array::toStdVector()
    {
        std::vector<double> vec = *this;
        return vec;
    }

    inline std::string Array::toString()
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

} // namespace Com
