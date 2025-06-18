#include "Modulator.h"
#include "SymbolMapper.h"

Com::Array<std::complex<double>> Com::Modulator::pskmod(const Array<double> &data, size_t M, double phase_offset, bool gray_coded)
{
    Array<std::complex<double>> modulated;
    modulated.reserve(data.size());

    double phase_inc = 2 * M_PI / M;

    for (double sym : data)
    {
        size_t idx = static_cast<size_t>(sym);
        if (idx >= M)
            throw std::invalid_argument("Symbol out of range");

        size_t mapped = gray_coded ? SymbolMapper::gray_encode(idx) : idx;
        double phase = mapped * phase_inc + phase_offset;
        modulated.push_back(std::polar(1.0, phase));
    }
    return modulated;
}

Com::Array<double> Com::Modulator::pskdemod(const Array<std::complex<double>> &modulated_signal, size_t M, bool gray_coded)
{
    Array<double> demodulated;
    double phase_inc = 2 * M_PI / M;
    double offset = phase_inc / 2.0;

    for (auto &z : modulated_signal)
    {
        double angle = std::arg(z);
        if (angle < 0)
            angle += 2 * M_PI;

        size_t idx = static_cast<size_t>(std::floor((angle + offset) / phase_inc)) % M;
        size_t decoded = gray_coded ? SymbolMapper::gray_decode(idx) : idx;
        demodulated.push_back(static_cast<double>(decoded));
    }
    return demodulated;
}

Com::Array<std::complex<double>> Com::Modulator::qammod(const Array<double> &data, size_t M, bool gray_coded)
{
    Array<std::complex<double>> modulated;
    modulated.reserve(data.size());

    std::vector<std::complex<double>> constellation;

    if (M == 32 || M == 128)
    {
        size_t k = (M == 32) ? 6 : 12; // Grid size: 6x6 for 32-QAM, 12x12 for 128-QAM
        int start = -(k - 1);          // Start value: -5 for M=32, -11 for M=128

        for (size_t i = 0; i < k; ++i)
        {
            for (size_t j = 0; j < k; ++j)
            {
                // Calculate coordinates (odd integers)
                double x = start + 2.0 * i;
                double y = start + 2.0 * j;

                // Remove corner points
                if (M == 32)
                {
                    // Remove 4 corners: (start,start), (start,-start), (-start,start), (-start,-start)
                    if ((i == 0 && j == 0) ||
                        (i == 0 && j == k - 1) ||
                        (i == k - 1 && j == 0) ||
                        (i == k - 1 && j == k - 1))
                    {
                        continue;
                    }
                }
                else if (M == 128)
                {
                    // Remove 16 corners: 2x2 blocks at each grid corner
                    if ((i < 2 && j < 2) ||      // Bottom-left block
                        (i < 2 && j >= k - 2) || // Top-left block
                        (i >= k - 2 && j < 2) || // Bottom-right block
                        (i >= k - 2 && j >= k - 2))
                    { // Top-right block
                        continue;
                    }
                }
                constellation.emplace_back(x, y);
            }
        }
    }
    else if (M == 8)
    {

        for (int y = -1; y <= 1; y += 2)
            for (int x = -3; x <= 3; x += 2)
                constellation.emplace_back(static_cast<double>(x), static_cast<double>(y));
    }
    else
    {
        size_t k = static_cast<size_t>(std::sqrt(M));
        if (k * k != M)
            throw std::invalid_argument("QAM constellation size must be a perfect square (except 8, 32, 128).");

        int half = static_cast<int>(k) / 2;
        for (int y = -half; y < static_cast<int>(k) - half; ++y)
            for (int x = -half; x < static_cast<int>(k) - half; ++x)
                constellation.emplace_back(2.0 * x, 2.0 * y);
    }

    // Normalize average energy to 1
    double energy = 0.0;
    for (const auto &pt : constellation)
        energy += std::norm(pt);
    double norm_factor = std::sqrt(energy / M);
    for (auto &pt : constellation)
        pt /= norm_factor;

    for (double symbol : data)
    {
        size_t index = static_cast<size_t>(symbol);
        if (index >= M)
            throw std::invalid_argument("Symbol index exceeds constellation size");

        size_t mapped = gray_coded ? SymbolMapper::gray_encode(index) : index;
        modulated.push_back(constellation[mapped]);
    }

    return modulated;
}

Com::Array<double> Com::Modulator::qamdemod(const Array<std::complex<double>> &modulated_signal, size_t M, bool gray_coded)
{
    Array<double> demodulated;
    std::vector<std::complex<double>> constellation;

    if (M == 32 || M == 128)
    {
        size_t k = (M == 32) ? 6 : 12; // Grid size: 6x6 for 32-QAM, 12x12 for 128-QAM
        int start = -(k - 1);          // Start value: -5 for M=32, -11 for M=128

        for (size_t i = 0; i < k; ++i)
        {
            for (size_t j = 0; j < k; ++j)
            {
                // Calculate coordinates (odd integers)
                double x = start + 2.0 * i;
                double y = start + 2.0 * j;

                // Remove corner points
                if (M == 32)
                {
                    // Remove 4 corners: (start,start), (start,-start), (-start,start), (-start,-start)
                    if ((i == 0 && j == 0) ||
                        (i == 0 && j == k - 1) ||
                        (i == k - 1 && j == 0) ||
                        (i == k - 1 && j == k - 1))
                    {
                        continue;
                    }
                }
                else if (M == 128)
                {
                    // Remove 16 corners: 2x2 blocks at each grid corner
                    if ((i < 2 && j < 2) ||      // Bottom-left block
                        (i < 2 && j >= k - 2) || // Top-left block
                        (i >= k - 2 && j < 2) || // Bottom-right block
                        (i >= k - 2 && j >= k - 2))
                    { // Top-right block
                        continue;
                    }
                }
                constellation.emplace_back(x, y);
            }
        }
    }
    else if (M == 8)
    {
        for (int y = -1; y <= 1; y += 2)
            for (int x = -3; x <= 3; x += 2)
                constellation.emplace_back(static_cast<double>(x), static_cast<double>(y));
    }
    else
    {
        size_t k = static_cast<size_t>(std::sqrt(M));
        if (k * k != M)
            throw std::invalid_argument("QAM constellation size must be a perfect square (except 8, 32, 128).");

        int half = static_cast<int>(k) / 2;
        for (int y = -half; y < static_cast<int>(k) - half; ++y)
            for (int x = -half; x < static_cast<int>(k) - half; ++x)
                constellation.emplace_back(2.0 * x, 2.0 * y);
    }

    // Normalize to unit average power
    double energy = 0.0;
    for (const auto &pt : constellation)
        energy += std::norm(pt);
    double norm_factor = std::sqrt(energy / M);
    for (auto &pt : constellation)
        pt /= norm_factor;

    for (const auto &z : modulated_signal)
    {
        size_t closest_idx = 0;
        double min_dist = std::norm(z - constellation[0]);
        for (size_t i = 1; i < constellation.size(); ++i)
        {
            double dist = std::norm(z - constellation[i]);
            if (dist < min_dist)
            {
                min_dist = dist;
                closest_idx = i;
            }
        }

        size_t decoded = gray_coded ? SymbolMapper::gray_decode(closest_idx) : closest_idx;
        demodulated.push_back(static_cast<double>(decoded));
    }

    return demodulated;
}

Com::Array<std::complex<double>> Com::Modulator::fskmod(const Array<double> &data, std::size_t M, double freq_sep, double nsamp, double Fs)
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

Com::Array<double> Com::Modulator::fskdemod(const Array<std::complex<double>> &modulated_signal, std::size_t M, double freq_sep, double nsamp, double Fs)
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

Com::Array<std::complex<double>> Com::Modulator::apskmod(const Array<double> &data, std::size_t M, const Array<double> &radii, const Array<std::complex<double>> &lut)
{
    if (lut.size() != M)
        throw std::invalid_argument("LUT must match constellation size");

    Array<std::complex<double>> modulated;
    modulated.reserve(data.size());

    for (double sym : data)
    {
        size_t idx = static_cast<size_t>(sym);
        if (idx >= M)
            throw std::invalid_argument("Symbol exceeds constellation size");
        modulated.push_back(lut[idx]);
    }

    return modulated;
}

Com::Array<double> Com::Modulator::apskdemod(const Array<std::complex<double>> &modulated_signal, std::size_t M, const Array<std::complex<double>> &lut)
{
    if (lut.size() != M)
        throw std::invalid_argument("LUT must match constellation size");

    Array<double> demodulated;
    demodulated.reserve(modulated_signal.size());

    for (auto &point : modulated_signal)
    {
        double min_dist = std::numeric_limits<double>::max();
        size_t best_idx = 0;

        for (size_t i = 0; i < M; ++i)
        {
            double dist = std::norm(point - lut[i]);
            if (dist < min_dist)
            {
                min_dist = dist;
                best_idx = i;
            }
        }
        demodulated.push_back(static_cast<double>(best_idx));
    }

    return demodulated;
}

Com::Array<std::complex<double>> Com::Modulator::dpskmod(const Array<double> &data, std::size_t M, double phase_offset)
{
    return Array<std::complex<double>>();
}

Com::Array<double> Com::Modulator::dpskdemod(const Array<std::complex<double>> &modulated_signal, std::size_t M)
{
    return Array<double>();
}

Com::Array<std::complex<double>> Com::Modulator::pammod(const Array<double> &data, std::size_t M, double phase_offset)
{
    Array<std::complex<double>> modulated_signal;
    double max_amplitude = (M - 1) / 2.0;
    double step_size = max_amplitude * 2 / (M - 1);

    for (double value : data)
    {
        double amplitude = (value * step_size) + max_amplitude;
        std::complex<double> modulated_point(amplitude, phase_offset);
        modulated_signal.push_back(modulated_point);
    }

    return modulated_signal;
}

Com::Array<double> Com::Modulator::pamdemod(const Array<std::complex<double>> &modulated_signal, std::size_t M)
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

Com::Array<std::complex<double>> Com::Modulator::pwmmod(const Array<double> &data, std::size_t M, double phase_offset)
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

Com::Array<double> Com::Modulator::pwmdemod(const Array<std::complex<double>> &modulated_signal, std::size_t M)
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

Com::Array<std::complex<double>> Com::Modulator::ppmmod(const Array<double> &data, std::size_t M, double phase_offset)
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

Com::Array<double> Com::Modulator::ppmdemod(const Array<std::complex<double>> &modulated_signal, std::size_t M)
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
