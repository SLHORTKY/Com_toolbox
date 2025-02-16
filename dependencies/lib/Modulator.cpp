#include "Modulator.h"

Com::Array<std::complex<double>> Com::Modulator::pskmod(const Array<double> &data, size_t M, double phase_offset)
{
    Array<std::complex<double>> vec;
    vec.reserve(data.size());
    double phase_increment = 2 * M_PI / M;

    for (size_t i = 0; i < data.size(); i++)
    {
        if (data[i] < 0 || data[i] >= M)
        {
            throw std::invalid_argument("Data values must be between 0 and M-1 for PSK modulation.");
        }

        double phase = phase_increment * data[i] + phase_offset;
        vec.push_back(std::polar(1.0, phase));
    }
    return vec;
}

Com::Array<double> Com::Modulator::pskdemod(const Com::Array<std::complex<double>> &modulated_signal, size_t M)
{
    Array<double> demodulated_data;
    demodulated_data.reserve(modulated_signal.size());

    double phase_increment = 2 * M_PI / M;
    double reference_offset = phase_increment / 2.0;  

    for (const auto &sample : modulated_signal)
    {
        double phase = std::arg(sample);  
       
        size_t symbol = static_cast<size_t>(std::floor((phase + reference_offset) / phase_increment)) % M;

        demodulated_data.push_back(static_cast<double>(symbol));
    }

    return demodulated_data;
}

Com::Array<std::complex<double>> Com::Modulator::qammod(const Com::Array<double> &data, size_t M, double phase_offset)
{
    if ((M & (M - 1)) != 0)
    {
        throw std::invalid_argument("M must be a power of 2");
    }

    int sqrtM = static_cast<int>(std::sqrt(M));
    double norm_factor = std::sqrt(2.0 / (3.0 * (M - 1))); 

    std::unordered_map<double, std::complex<double>> qam_map;
    Com::Array<std::complex<double>> symbols;

    for (int i = 0; i < sqrtM; ++i)
    {
        for (int j = 0; j < sqrtM; ++j)
        {
            double I = (2 * i - (sqrtM - 1));
            double Q = (2 * j - (sqrtM - 1));
            std::complex<double> qamSymbol = std::complex<double>(I, Q) * norm_factor;  
            qamSymbol *= std::polar(1.0, phase_offset);  
            qam_map[i * sqrtM + j] = qamSymbol;
        }
    }

    for (double d : data)
    {
        symbols.push_back(qam_map[d]);
    }
    return symbols;
}

Com::Array<double> Com::Modulator::qamdemod(const Com::Array<std::complex<double>> &modulated_signal, size_t M)
{
    if ((M & (M - 1)) != 0 )
    {
        throw std::invalid_argument("M must be a power of 2");
    }

    int sqrtM = static_cast<int>(std::sqrt(M));
    double norm_factor = std::sqrt(2.0 / (3.0 * (M - 1)));  

    std::vector<std::complex<double>> qamSymbols;
    for (int i = 0; i < sqrtM; ++i)
    {
        for (int j = 0; j < sqrtM; ++j)
        {
            double I = (2 * i - (sqrtM - 1));
            double Q = (2 * j - (sqrtM - 1));
            qamSymbols.push_back(std::complex<double>(I, Q) * norm_factor);  // ✅ Normalize
        }
    }

    Com::Array<double> symbols;
    for (const auto &received_symbol : modulated_signal)
    {
        double min_distance = std::numeric_limits<double>::max();
        int closest_symbol = 0;

        for (size_t k = 0; k < qamSymbols.size(); ++k)
        {
            double distance = std::norm(received_symbol - qamSymbols[k]);  // ✅ Euclidean distance squared
            if (distance < min_distance)
            {
                min_distance = distance;
                closest_symbol = static_cast<int>(k);
            }
        }
        symbols.push_back(closest_symbol);
    }

    return symbols;
}

Com::Array<std::complex<double>> Com::Modulator::fskmod(const Array<double> &data, size_t M, double phase_offset)
{
    return Array<std::complex<double>>();
}

Com::Array<double> Com::Modulator::fskdemod(const Array<std::complex<double>> &modulated_signal, size_t M)
{
    return Array<double>();
}

Com::Array<std::complex<double>> Com::Modulator::apskmod(const Array<double> &data, size_t M, double phase_offset)
{
    return Array<std::complex<double>>();
}

Com::Array<double> Com::Modulator::apskdemod(const Array<std::complex<double>> &modulated_signal, size_t M)
{
    return Array<double>();
}

Com::Array<std::complex<double>> Com::Modulator::dpskmod(const Array<double> &data, size_t M, double phase_offset)
{
    return Array<std::complex<double>>();
}

Com::Array<double> Com::Modulator::dpskdemod(const Array<std::complex<double>> &modulated_signal, size_t M)
{
    return Array<double>();
}

Com::Array<std::complex<double>> Com::Modulator::pammod(const Array<double> &data, size_t M, double phase_offset)
{
    return Array<std::complex<double>>();
}

Com::Array<double> Com::Modulator::pamdemod(const Array<std::complex<double>> &modulated_signal, size_t M)
{
    return Array<double>();
}

Com::Array<std::complex<double>> Com::Modulator::pwmmod(const Array<double> &data, size_t M, double phase_offset)
{
    return Array<std::complex<double>>();
}

Com::Array<double> Com::Modulator::pwmdemod(const Array<std::complex<double>> &modulated_signal, size_t M)
{
    return Array<double>();
}

Com::Array<std::complex<double>> Com::Modulator::ppmmod(const Array<double> &data, size_t M, double phase_offset)
{
    return Array<std::complex<double>>();
}

Com::Array<double> Com::Modulator::ppmdemod(const Array<std::complex<double>> &modulated_signal, size_t M)
{
    return Array<double>();
}
