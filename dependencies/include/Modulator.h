#pragma once
#include "Array.h"

namespace Com
{
    class Modulator
    {
    public:
        // Phase Shift Keying (PSK) and Demodulation (PSK)
        static Array<std::complex<double>> pskmod(const Array<double>& data, std::size_t M, double phase_offset, bool gray_coded = true);
        static Array<double> pskdemod(const Array<std::complex<double>>& modulated_signal, std::size_t M, bool gray_coded = true);

        // Quadrature Amplitude Modulation (QAM) and Demodulation (QAM)
        static Array<std::complex<double>> qammod(const Array<double>& data, std::size_t M, bool gray_coded = true);
        static Array<double> qamdemod(const Array<std::complex<double>>& modulated_signal, std::size_t M, bool gray_coded = true);

        // Frequency Shift Keying Modulation (FSK) and Demodulation (FSK)
        static Array<std::complex<double>> fskmod(const Array<double>& data, std::size_t M, double freq_sep, double nsamp, double Fs);
        static Array<double> fskdemod(const Array<std::complex<double>>& modulated_signal, std::size_t M, double freq_sep, double nsamp, double Fs);

        // Amplitude Phase Shift Keying Modulation (APSK) and Demodulation (APSK)
        static Array<std::complex<double>> apskmod(const Array<double>& data, std::size_t M, const Array<double>& radii, const Array<std::complex<double>>& lut);
        static Array<double> apskdemod(const Array<std::complex<double>>& modulated_signal, std::size_t M, const Array<std::complex<double>>& lut);

        // Differential Phase Shift Keying Modulation (DPSK) and Demodulation (DPSK)
        static Array<std::complex<double>> dpskmod(const Array<double>& data, std::size_t M, double phase_offset);
        static Array<double> dpskdemod(const Array<std::complex<double>>& modulated_signal, std::size_t M);

        // Pulse Amplitude Modulation (PAM) and Demodulation (PAM)
        static Array<std::complex<double>> pammod(const Array<double>& data, std::size_t M, double phase_offset);
        static Array<double> pamdemod(const Array<std::complex<double>>& modulated_signal, std::size_t M);

        // Pulse Width Modulation (PWM) and Demodulation (PWM)
        static Array<std::complex<double>> pwmmod(const Array<double>& data, std::size_t M, double phase_offset);
        static Array<double> pwmdemod(const Array<std::complex<double>>& modulated_signal, std::size_t M);

        // Pulse Position Modulation (PPM) and Demodulation (PPM)
        static Array<std::complex<double>> ppmmod(const Array<double>& data, std::size_t M, double phase_offset);
        static Array<double> ppmdemod(const Array<std::complex<double>>& modulated_signal, std::size_t M);

    };

} // namespace Com
