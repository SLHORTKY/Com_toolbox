#include "Filter.h"

using namespace Com;

Array<double> Com::Filter::generateSquareFilter(int length, double amplitude)
{
    return Array<double>::ones(length) * amplitude;
}
// RC Filter
Array<double> Com::Filter::generateRaisedCosineFilter(double beta, double sps, double span)
{
    if (beta < 0.0 || beta > 1.0)
    {
        throw std::invalid_argument("Beta must be in the range [0, 1].");
    }
    if (sps <= 0 || span <= 0)
    {
        throw std::invalid_argument("Samples per symbol (sps) and span must be positive.");
    }

    size_t num_taps = static_cast<size_t>(span * sps) + 1;
    
    Array<double> h = Array<double>::ones(num_taps) * 0.0;

    double T = 1.0;
    double Ts = T / sps;              
    int t_mid = (num_taps - 1) / 2.0; 

    for (int i = 0; i < num_taps; ++i)
    {
        double t = (i - t_mid) * Ts; // Time index

        if (t == 0.0)
        {
            h[i] = 1.0;
        }
        else if (beta != 0.0 && std::abs(t) == (T / (2 * beta)))
        {
            h[i] = (beta / std::sqrt(2)) * ((M_PI / 4) * std::sin(M_PI / (2 * beta)));
        }
        else
        {
            double sinc_val = std::sin(M_PI * t / T) / (M_PI * t / T);
            double cos_val = std::cos(M_PI * beta * t / T);
            double denom = 1 - (4 * beta * beta * t * t) / (T * T);
            h[i] = sinc_val * cos_val / denom;
        }
    }

    double sum = h.apply(SignalMath::sum);
    Array<double> normalized_h = h.apply(SignalMath::normalize, sum);
    return normalized_h;
}
// RRC Filter
Array<double> Com::Filter::generateRootRaisedCosineFilter(double beta, double sps, double span)
{
    return generateRaisedCosineFilter(beta, sps, span).apply(SignalMath::abs).apply(SignalMath::sqrt);
}
// Gaussian Filter
Array<double> Com::Filter::generateGaussianFilter(int length, double bandwidth)
{
    if (length <= 0 || bandwidth <= 0)
    {
        throw std::invalid_argument("Length and bandwidth must be positive.");
    }

    Array<double> h(length);
    double sigma = bandwidth / std::sqrt(2 * std::log(2)); // Convert bandwidth to standard deviation
    int mid = length / 2;
    double sum = 0.0;

    for (int i = 0; i < length; i++)
    {
        double x = i - mid;
        h[i] = std::exp(-0.5 * (x * x) / (sigma * sigma));
        sum += h[i];
    }

    Array<double> h_normalized = h.apply(SignalMath::normalize, sum);
    return h_normalized;
}
// Butterworth Filter
Array<double> Com::Filter::generateButterwortFilter(int n, double wc, std::string type)
{
    if (n <= 0 || wc <= 0 || wc >= 1)
    {
        throw std::invalid_argument("Order must be positive and cutoff frequency in (0,1).");
    }

    int length = 2 * n + 1; // Number of taps
    Array<double> h(length);
    int mid = length / 2;
    double sum = 0.0;

    for (int i = 0; i < length; i++)
    {
        double x = i - mid;
        double omega = wc * M_PI;
        double butterworth = 1.0 / std::sqrt(1 + std::pow(x / omega, 2 * n));

        if (type == "low")
        {
            h[i] = butterworth;
        }
        else if (type == "high")
        {
            h[i] = 1 - butterworth;
        }
        else if (type == "stop")
        {
            h[i] = 1 - 2 * butterworth;
        }
        else
        {
            throw std::invalid_argument("Invalid filter type. Use 'low', 'high', or 'stop'.");
        }

        sum += h[i];
    }

    Array<double> h_normalized = h.apply(SignalMath::normalize, sum);
    return h_normalized;
}
// Chebyshev Filter
Array<double> Com::Filter::generateChebyshevFilter(int n, double ripple, double wc)
{
    if (n <= 0 || ripple <= 0 || wc <= 0 || wc >= 1)
    {
        throw std::invalid_argument("Invalid filter parameters.");
    }

    int length = 2 * n + 1;
    Array<double> h(length);
    int mid = length / 2;
    double sum = 0.0;
    double epsilon = std::sqrt(std::pow(10, ripple / 10.0) - 1);
    double omega_c = wc * M_PI;

    for (int i = 0; i < length; i++)
    {
        double x = i - mid;
        double Tn = std::cosh(n * std::acosh(x / omega_c)); // Chebyshev polynomial
        h[i] = 1.0 / std::sqrt(1 + epsilon * epsilon * Tn * Tn);
        sum += h[i];
    }

    Array<double> h_normalized = h.apply(SignalMath::normalize, sum);
    return h_normalized;
}
// Matched Filter
Array<double> Com::Filter::generateMatchedFilter(const Array<double> &inputSignal)
{
    int len = inputSignal.size();
    if (len == 0)
        return Array<double>();

    Array<double> h(len);
    for (int i = 0; i < len; i++)
    {
        h[i] = inputSignal[len - 1 - i]; // Time-reversed signal
    }
    return h;
}
// FIR Filter
Array<double> Com::Filter::generateFIRFilter(const Array<double> &coefficients, size_t numSamples)
{
    if (coefficients.empty())
    {
        throw std::invalid_argument("FIR filter coefficients cannot be empty.");
    }

    Array<double> impulseResponse = Array<double>::ones(numSamples) * 0.0;

    for (size_t n = 0; n < coefficients.size(); ++n)
    {
        impulseResponse[n] = coefficients[n];
    }

    return impulseResponse;
}
// IIR Filter
Array<double> Com::Filter::generateIIRFilter(const Array<double> &aCoeffs, const Array<double> &bCoeffs, size_t numSamples)
{
    if (aCoeffs.empty() || bCoeffs.empty())
    {
        throw std::invalid_argument("IIR filter coefficients cannot be empty.");
    }

    Array<double> impulseResponse = Array<double>::ones(numSamples) * 0.0;

    Array<double> x = Array<double>::ones(numSamples) * 0.0;
    x[0] = 1.0;

    Array<double> y = Array<double>::ones(numSamples) * 0.0;

    for (size_t n = 0; n < numSamples; ++n)
    {
        for (size_t i = 0; i < bCoeffs.size(); ++i)
        {
            if (n >= i)
                y[n] += bCoeffs[i] * x[n - i];
        }

        for (size_t i = 1; i < aCoeffs.size(); ++i)
        {
            if (n >= i)
                y[n] -= aCoeffs[i] * y[n - i];
        }

        impulseResponse[n] = y[n];
    }

    return impulseResponse;
}
// LMS Adaptive Filter
Array<double> Com::Filter::applyLMSFilter(const Array<double> &inputSignal, const Array<double> &desiredSignal, double stepSize)
{
    if (inputSignal.size() != desiredSignal.size())
    {
        throw std::invalid_argument("Input and desired signals must have the same length.");
    }

    int len = inputSignal.size();
    Array<double> weights = Array<double>::ones((len)) * 0.0;
    Array<double> output = Array<double>::ones((len)) * 0.0;

    for (int i = 0; i < len; i++)
    {
        double error = desiredSignal[i] - output[i];   
        weights[i] += stepSize * error * inputSignal[i]; 
        output[i] = weights[i] * inputSignal[i];        
    }

    return output;
}
// Comb Filter (Echo Effect)
Array<double> Com::Filter::generateCombFilter(int delay, double gain)
{
    if (delay <= 0 || gain <= 0)
    {
        throw std::invalid_argument("Delay and gain must be positive.");
    }

    Array<double> h = Array<double>::ones((delay + 1)) * 0.0;
    h[0] = 1.0;
    h[delay] = gain; // Delayed echo

    return h;
}
