#pragma once
#include "Array.h"
#include "Filter.h"
#include "SignalMath.h"

namespace Com
{
    template <typename NUMERIC>
    class ComSignal
    {
    private:
        Array<NUMERIC> signal;
        Array<double> domain;

    public:
        ComSignal(Array<NUMERIC> signal, Array<double> domain);
        ComSignal(const std::vector<NUMERIC>& vec) : signal(vec) {}

        ~ComSignal() = default;

        ComSignal convolution(const ComSignal &other) const;
        ComSignal crossCorrelation(const ComSignal &other) const;
        ComSignal autoCorrelation() const;
        ComSignal fourierTransform() const;
        ComSignal applyFilter(const Filter &filter) const;
        ComSignal upSample(size_t factor) const;
        ComSignal downSample(size_t factor) const;

        size_t size() const;
        Array<NUMERIC> Signal() const;
        Array<double> Domain() const;

        void setSignal(Array<NUMERIC> value);
        void setDomain(Array<double> value);
        
        double energy() const;
        double power() const;
    };

    template <typename NUMERIC>
    inline ComSignal<NUMERIC>::ComSignal(Array<NUMERIC> signal, Array<double> domain)
    {
        this->signal = signal;
        this->domain = domain;
    }
    template <typename NUMERIC>
    inline ComSignal<NUMERIC> ComSignal<NUMERIC>::convolution(const ComSignal &other) const
    {

        size_t len1 = this->size();
        size_t len2 = other.size();

        std::vector<NUMERIC> result(len1 + len2 - 1, 0);

        for (size_t i = 0; i < len1; ++i)
        {
            for (size_t j = 0; j < len2; ++j)
            {
                result[i + j] += this->signal[i] * other.signal[j];
            }
        }

        return ComSignal<NUMERIC>(result, Array<double>::linespace(0,result.size(),result.size()));
    }
    template <typename NUMERIC>
    inline ComSignal<NUMERIC> ComSignal<NUMERIC>::crossCorrelation(const ComSignal &other) const
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
                    result[m] += this->signal[n] * other.signal[m - n];
                }
            }
        }

        return ComSignal<NUMERIC>(result);
    }
    template <typename NUMERIC>
    inline ComSignal<NUMERIC> ComSignal<NUMERIC>::autoCorrelation() const
    {
        return this->crossCorrelation(*this);
    }
    template <typename NUMERIC>
    inline ComSignal<NUMERIC> ComSignal<NUMERIC>::fourierTransform() const
    {
        return ComSignal();
    }
    template <typename NUMERIC>
    inline ComSignal<NUMERIC> ComSignal<NUMERIC>::applyFilter(const Filter &filter) const
    {
        return this->convolution(filter);
    }
    template <typename NUMERIC>
    inline ComSignal<NUMERIC> ComSignal<NUMERIC>::upSample(size_t factor) const
    {
        if (factor < 1)
        {
            throw std::invalid_argument("Upsampling factor must be at least 1.");
        }

        std::vector<NUMERIC> result;
        result.reserve(this->size() * factor);

        for (size_t i = 0; i < this->size(); ++i)
        {
            result.push_back(this->data()[i]); // Original sample
            for (size_t j = 1; j < factor; ++j)
            {
                result.push_back(0); // Insert zeros
            }
        }

        return ComSignal<NUMERIC>(result);
    }
    template <typename NUMERIC>
    inline ComSignal<NUMERIC> ComSignal<NUMERIC>::downSample(size_t factor) const
    {
        if (factor < 1)
        {
            throw std::invalid_argument("Downsampling factor must be at least 1.");
        }

        std::vector<NUMERIC> result;
        result.reserve(this->size() / factor);

        for (size_t i = 0; i < this->size(); i += factor)
        {
            result.push_back(this->data()[i]); // Keep every `factor`-th sample
        }

        return ComSignal<NUMERIC>(result);
    }
    template <typename NUMERIC>
    inline size_t ComSignal<NUMERIC>::size() const
    {
        return this->signal.size();
    }
    template <typename NUMERIC>
    inline Array<NUMERIC> ComSignal<NUMERIC>::Signal() const
    {
        return this->signal;
    }
    template <typename NUMERIC>
    inline Array<double> ComSignal<NUMERIC>::Domain() const
    {
        return this->domain;
    }
    template <typename NUMERIC>
    inline void ComSignal<NUMERIC>::setSignal(Array<NUMERIC> value)
    {
        this->signal = value;
    }

    template <typename NUMERIC>
    inline void ComSignal<NUMERIC>::setDomain(Array<double> value)
    {
        this->domain = value;
    }
    template <typename NUMERIC>
    inline double ComSignal<NUMERIC>::energy() const
    {
        return this->signal.apply(SignalMath::pow, 2.0).apply(SignalMath::sum);
    }
    template <typename NUMERIC>
    inline double ComSignal<NUMERIC>::power() const
    {
        return this->signal.apply(SignalMath::pow, 2.0).apply(SignalMath::mean);
    }

} // namespace Com
