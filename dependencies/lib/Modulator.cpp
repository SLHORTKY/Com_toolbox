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
}

Com::Array<double> Com::Modulator::qamdemod(const Com::Array<std::complex<double>> &modulated_signal, size_t M)
{
}

Com::Array<std::complex<double>> Com::Modulator::fskmod(const Array<double> &data, size_t M, double freq_sep, double nsamp, double Fs)
{
    Array<std::complex<double>> modulated_signal;
    size_t samples_per_symbol = static_cast<size_t>(nsamp);
    double two_pi = 2.0 * M_PI;

    for (double symbol : data)
    {
        double freq = (symbol - (M / 2)) * freq_sep; // Centered frequencies
        for (size_t n = 0; n < samples_per_symbol; ++n)
        {
            double t = static_cast<double>(n) / Fs;
            std::complex<double> sample = std::polar(1.0, two_pi * freq * t);
            modulated_signal.push_back(sample);
        }
    }
    return modulated_signal;
}

Com::Array<double> Com::Modulator::fskdemod(const Array<std::complex<double>> &modulated_signal, size_t M, double freq_sep, double nsamp, double Fs)
{
    Array<double> demodulated_data;
    size_t samples_per_symbol = static_cast<size_t>(nsamp);
    size_t num_symbols = modulated_signal.size() / samples_per_symbol;
    double two_pi = 2.0 * M_PI;

    for (size_t i = 0; i < num_symbols; ++i)
    {
        std::vector<double> power(M, 0.0);
        for (size_t m = 0; m < M; ++m)
        {
            double freq = (m - (M / 2)) * freq_sep;
            std::complex<double> sum(0.0, 0.0);
            for (size_t n = 0; n < samples_per_symbol; ++n)
            {
                double t = static_cast<double>(n) / Fs;
                std::complex<double> reference = std::polar(1.0, -two_pi * freq * t);
                sum += modulated_signal[i * samples_per_symbol + n] * reference;
            }
            power[m] = std::norm(sum);
        }
        demodulated_data.push_back(std::distance(power.begin(), std::max_element(power.begin(), power.end())));
    }
    return demodulated_data;
}

Com::Array<std::complex<double>> Com::Modulator::apskmod(const Array<double> &data, size_t M, const Array<double> &radii)
{
    return Array<std::complex<double>>();
}

Com::Array<double> Com::Modulator::apskdemod(const Array<std::complex<double>>& modulated_signal, size_t M , const Array<double> &radii)
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
    Array<std::complex<double>> modulated_signal;
    double max_amplitude = (M - 1) / 2.0;
    double step_size = max_amplitude * 2 / (M - 1);

    // Modulating each data symbol
    for (double value : data)
    {
        // Find the corresponding amplitude
        double amplitude = (value * step_size) + max_amplitude;
        // Create complex signal (real part is the amplitude, imaginary part can be phase offset)
        std::complex<double> modulated_point(amplitude, phase_offset);
        modulated_signal.push_back(modulated_point);
    }

    return modulated_signal;
}

Com::Array<double> Com::Modulator::pamdemod(const Array<std::complex<double>> &modulated_signal, size_t M)
{
    Array<double> demodulated_data;
    double max_amplitude = (M - 1) / 2.0;
    double step_size = max_amplitude * 2 / (M - 1);

    // Demodulating each signal point
    for (const auto &modulated_point : modulated_signal)
    {
        double received_amplitude = modulated_point.real();
        double symbol = std::round((received_amplitude - max_amplitude) / step_size);
        symbol = std::clamp(symbol, 0.0, static_cast<double>(M - 1));
        demodulated_data.push_back(symbol);
    }

    return demodulated_data;
}

Com::Array<std::complex<double>> Com::Modulator::pwmmod(const Array<double> &data, size_t M, double phase_offset)
{
    Array<std::complex<double>> modulated_signal;
    double max_width = (M - 1);

    // Modulating each data symbol
    for (double value : data)
    {
        // Map to pulse width
        double width = (value * max_width);
        std::complex<double> modulated_point(width, phase_offset);
        modulated_signal.push_back(modulated_point);
    }

    return modulated_signal;
}

Com::Array<double> Com::Modulator::pwmdemod(const Array<std::complex<double>> &modulated_signal, size_t M)
{
    Array<double> demodulated_data;
    double max_width = (M - 1);

    // Demodulating each signal point
    for (const auto &modulated_point : modulated_signal)
    {
        double pulse_width = modulated_point.real();
        double symbol = pulse_width / max_width;
        demodulated_data.push_back(symbol);
    }

    return demodulated_data;
}

Com::Array<std::complex<double>> Com::Modulator::ppmmod(const Array<double> &data, size_t M, double phase_offset)
{
    Array<std::complex<double>> modulated_signal;
    double max_position = 1.0; // Maximum position (scaled)
    double step_size = max_position / M;

    // Modulating each data symbol
    for (double value : data)
    {
        // Map to pulse position
        double position = value * step_size;
        std::complex<double> modulated_point(position, phase_offset);
        modulated_signal.push_back(modulated_point);
    }

    return modulated_signal;
}

Com::Array<double> Com::Modulator::ppmdemod(const Array<std::complex<double>> &modulated_signal, size_t M)
{
    Array<double> demodulated_data;
    double max_position = 1.0;
    double step_size = max_position / M;

    // Demodulating each signal point
    for (const auto &modulated_point : modulated_signal)
    {
        double position = modulated_point.real();
        double symbol = position / step_size;
        demodulated_data.push_back(symbol);
    }

    return demodulated_data;
}
