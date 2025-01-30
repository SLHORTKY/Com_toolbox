#pragma once
#include <iostream>
#include "Array.h"
#include <random>
#include <vector>
#include <functional>
#include <unordered_map>

namespace Com
{
    class Noise
    {
    private:
        std::random_device rd;
        std::mt19937 gen;

    public:
        Noise() : gen(rd()) {}
        ~Noise() = default;

        Array<double> generateGaussian(double mean, double stddev, std::size_t n);
        Array<double> generateUniform(double start, double end, std::size_t n);
        Array<double> generatePoisson(double mean, std::size_t n);
        Array<double> generateBernoulli(double p, std::size_t n);
        Array<double> generateBinomial(int trials, double p, std::size_t n);
        Array<double> generateExponential(double lambda, std::size_t n);
        Array<double> generateWeibull(double shape, double scale, std::size_t n);
        Array<double> generateLognormal(double mean, double stddev, std::size_t n);
        Array<double> generateChiSquared(double df, std::size_t n);
        Array<double> generateStudent_t(double df, std::size_t n);
        Array<double> generateDiscrete(const std::vector<double> &probabilities, std::size_t n);
        Array<double> generateCauchy(double alpha, double beta, std::size_t n);
        Array<double> generateGamma(double alpha, double beta, std::size_t n);
    };

    // Function to generate n random numbers in the range [start, end]
    inline Array<double> Noise::generateUniform(double start, double end, std::size_t n)
    {
        std::uniform_real_distribution<double> dis(start, end);
        std::vector<double> result(n);
        for (std::size_t i = 0; i < n; ++i)
        {
            result[i] = dis(gen);
        }
        return Array<double>(result);
    }
    // Gaussian (Normal) distribution
    inline Array<double> Noise::generateGaussian(double mean, double stddev, std::size_t n)
    {
        std::normal_distribution<double> dis(mean, stddev);
        std::vector<double> result(n);
        for (std::size_t i = 0; i < n; ++i)
        {
            result[i] = dis(gen);
        }
        return Array<double>(result);
    }
    // Poisson distribution
    inline Array<double> Noise::generatePoisson(double mean, std::size_t n)
    {
        std::poisson_distribution<int> dis(mean);
        std::vector<double> result(n);
        for (std::size_t i = 0; i < n; ++i)
        {
            result[i] = dis(gen);
        }
        return Array<double>(result);
    }
    // Bernoulli distribution
    inline Array<double> Noise::generateBernoulli(double p, std::size_t n)
    {
        std::bernoulli_distribution dis(p);
        std::vector<double> result(n);
        for (std::size_t i = 0; i < n; ++i)
        {
            result[i] = dis(gen) ? 1.0 : 0.0;
        }
        return Array<double>(result);
    }
    // Binomial distribution
    inline Array<double> Noise::generateBinomial(int trials, double p, std::size_t n)
    {
        std::binomial_distribution<int> dis(trials, p);
        std::vector<double> result(n);
        for (std::size_t i = 0; i < n; ++i)
        {
            result[i] = dis(gen);
        }
        return Array<double>(result);
    }
    // Exponential distribution
    inline Array<double> Noise::generateExponential(double lambda, std::size_t n)
    {
        std::exponential_distribution<double> dis(lambda);
        std::vector<double> result(n);
        for (std::size_t i = 0; i < n; ++i)
        {
            result[i] = dis(gen);
        }
        return Array<double>(result);
    }
    // Gamma distribution
    inline Array<double> Noise::generateGamma(double alpha, double beta, std::size_t n)
    {
        std::gamma_distribution<double> dis(alpha, beta);
        std::vector<double> result(n);
        for (std::size_t i = 0; i < n; ++i)
        {
            result[i] = dis(gen);
        }
        return Array<double>(result);
    }
    // Weibull distribution
    inline Array<double> Noise::generateWeibull(double shape, double scale, std::size_t n)
    {
        std::weibull_distribution<double> dis(shape, scale);
        std::vector<double> result(n);
        for (std::size_t i = 0; i < n; ++i)
        {
            result[i] = dis(gen);
        }
        return Array<double>(result);
    }
    // Lognormal distribution
    inline Array<double> Noise::generateLognormal(double mean, double stddev, std::size_t n)
    {
        std::lognormal_distribution<double> dis(mean, stddev);
        std::vector<double> result(n);
        for (std::size_t i = 0; i < n; ++i)
        {
            result[i] = dis(gen);
        }
        return Array<double>(result);
    }
    // Chi-squared distribution
    inline Array<double> Noise::generateChiSquared(double df, std::size_t n)
    {
        std::chi_squared_distribution<double> dis(df);
        std::vector<double> result(n);
        for (std::size_t i = 0; i < n; ++i)
        {
            result[i] = dis(gen);
        }
        return Array<double>(result);
    }
    // Student's t-distribution
    inline Array<double> Noise::generateStudent_t(double df, std::size_t n)
    {
        std::student_t_distribution<double> dis(df);
        std::vector<double> result(n);
        for (std::size_t i = 0; i < n; ++i)
        {
            result[i] = dis(gen);
        }
        return Array<double>(result);
    }
    // Beta distribution
    inline Array<double> Noise::generateCauchy(double alpha, double beta, std::size_t n)
    {
        std::cauchy_distribution<double> dis(alpha, beta);
        std::vector<double> result(n);
        for (std::size_t i = 0; i < n; ++i)
        {
            result[i] = dis(gen);
        }
        return Array<double>(result);
    }
    // Discrete distribution
    inline Array<double> Noise::generateDiscrete(const std::vector<double> &probabilities, std::size_t n)
    {
        // Create a local copy of the probabilities vector to modify (avoid modifying the original)
        std::vector<double> modifiedProbabilities = probabilities;

        // Check if probabilities sum to 1.0 (optional)
        double sum = std::accumulate(modifiedProbabilities.begin(), modifiedProbabilities.end(), 0.0);
        if (sum != 1.0)
        {
            // Normalize the probabilities if they do not sum to 1
            for (auto &p : modifiedProbabilities)
            {
                p /= sum;
            }
        }
        // Create the discrete distribution using the modified probabilities
        std::discrete_distribution<int> dis(modifiedProbabilities.begin(), modifiedProbabilities.end());

        // Create a result vector of size n
        std::vector<double> result(n);

        // Generate n samples based on the distribution
        for (std::size_t i = 0; i < n; ++i)
        {
            result[i] = static_cast<double>(dis(gen)); // Convert the result to double
        }

        // Return the result wrapped in Array<double>
        return Array<double>(result);
    }

} // namespace Com
