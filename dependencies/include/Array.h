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

namespace Com
{
    template <typename NUMERIC> // INT, //DOUBLE //SHORT ... std::complex
    class Array : public std::vector<NUMERIC>
    {
    public:
        Array(std::initializer_list<NUMERIC> init);
        Array(const std::vector<NUMERIC> &vec);
        Array(size_t size);
        Array() = default;

        static Array<double> arange(NUMERIC start, NUMERIC end, NUMERIC step);
        static Array<double> linespace(NUMERIC start, NUMERIC end, size_t num);

        virtual Array operator+(const Array &other) const;
        virtual Array operator-(const Array &other) const;
        virtual Array operator*(const Array &other) const;
        virtual Array operator/(const Array &other) const;

        virtual Array operator+(double other) const;
        virtual Array operator-(double other) const;
        virtual Array operator*(double other) const;
        virtual Array operator/(double other) const;

        Array slice(size_t start, size_t end, size_t step) const;
        Array operator()(size_t start, size_t end, size_t step) const;

        Array apply(std::function<NUMERIC(NUMERIC)> func) const;
        Array apply(NUMERIC (*func)(NUMERIC)) const;
        Array apply(NUMERIC (*func)(NUMERIC, NUMERIC), const Array &other) const;
        NUMERIC apply(NUMERIC (*func)(const std::vector<NUMERIC> &vec)) const;
        Array apply(NUMERIC (*func)(NUMERIC, std::unordered_map<std::string, NUMERIC> &params), std::unordered_map<std::string, NUMERIC> &params) const;

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
    };
    template <typename NUMERIC>
    inline Array<NUMERIC>::Array(std::initializer_list<NUMERIC> init)
        : std::vector<NUMERIC>(init) {}

    template <typename NUMERIC>
    inline Array<NUMERIC>::Array(const std::vector<NUMERIC> &vec)
        : std::vector<NUMERIC>(vec) {}

    template <typename NUMERIC>
    inline Array<NUMERIC>::Array(size_t size)
        : std::vector<NUMERIC>(size) {}

    template <typename NUMERIC>
    inline Array<double> Array<NUMERIC>::arange(NUMERIC start, NUMERIC end, NUMERIC step)
    {
        static_assert(std::is_same_v<NUMERIC, double> || std::is_same_v<NUMERIC, int>,
                      "arange is only supported for int and double types");

        if (step == 0)
            throw std::invalid_argument("Step cannot be zero");

        Array result;
        for (double value = start; (step > 0 ? value < end : value > end); value += step)
        {
            result.push_back(value);
        }
        return result;
    }
    template <typename NUMERIC>
    inline Array<double> Array<NUMERIC>::linespace(NUMERIC start, NUMERIC end, size_t num)
    {
        static_assert(std::is_same_v<NUMERIC, double> || std::is_same_v<NUMERIC, int>,
                      "arange is only supported for int and double types");

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
    inline Array<NUMERIC> Array<NUMERIC>::operator+(double other) const
    {
        Array result(this->size());
        for (size_t i = 0; i < this->size(); ++i)
        {
            result[i] = this->at(i) + other;
        }
        return result;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::operator-(double other) const
    {
        Array result(this->size());
        for (size_t i = 0; i < this->size(); ++i)
        {
            result[i] = this->at(i) - other;
        }
        return result;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::operator*(double other) const
    {
        Array result(this->size());
        for (size_t i = 0; i < this->size(); ++i)
        {
            result[i] = this->at(i) * other;
        }
        return result;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::operator/(double other) const
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
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::operator()(size_t start, size_t end, size_t step) const
    {
        return this->slice(start, end, step);
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::slice(size_t start, size_t end, size_t step) const
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
        Array vec;
        for (size_t i = 0; i < this->size(); i++)
        {
            vec.push_back(func(this->at(i)));
        }
        return vec;
    }
    template <typename NUMERIC>
    inline NUMERIC Array<NUMERIC>::apply(NUMERIC (*func)(const std::vector<NUMERIC> &vec)) const
    {
        return func(*this);
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::apply(NUMERIC (*func)(NUMERIC)) const
    {
        Array vec;
        for (size_t i = 0; i < this->size(); i++)
        {
            vec.push_back(func(this->at(i)));
        }
        return vec;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::apply(NUMERIC (*func)(NUMERIC, NUMERIC), const Array &other) const
    {
        Array vec;
        for (size_t i = 0; i < this->size(); i++)
        {
            vec.push_back(func(this->at(i), other.at(i)));
        }
        return vec;
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> Array<NUMERIC>::apply(NUMERIC (*func)(NUMERIC, std::unordered_map<std::string, NUMERIC> &params), std::unordered_map<std::string, NUMERIC> &params) const
    {
        std::vector<double> vec;
        for (size_t i = 0; i < this->size(); i++)
        {
            vec.push_back(func(this->at(i), params));
        }
        return vec;
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

} // namespace Com
