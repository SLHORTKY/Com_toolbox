#pragma once
#include "Array.h"
#include <iostream>

namespace Com
{
    class Modulator
    {
    public:
        static Array<std::complex<double>> pskmod(Array<double> data, int16_t M, double phase_offset);
        static Array<double> pskdemod(Array<std::complex<double>> modulated_signal, int16_t M);

        // Quadrature Amplitude Modulation (QAM) and Demodulation (QAM)
        static Array<std::complex<double>> qammod(Array<double> data, int16_t M, double phase_offset);
        static Array<double> qamdemod(Array<std::complex<double>> modulated_signal, int16_t M);

        // Frequency Shift Keying Modulation (FSK) and Demodulation (FSK)
        static Array<std::complex<double>> fskmod(Array<double> data, int16_t M, double phase_offset);
        static Array<double> fskdemod(Array<std::complex<double>> modulated_signal, int16_t M);

        // Amplitude Phase Shift Keying Modulation (APSK) and Demodulation (APSK)
        static Array<std::complex<double>> apskmod(Array<double> data, int16_t M, double phase_offset);
        static Array<double> apskdemod(Array<std::complex<double>> modulated_signal, int16_t M);

        // Differential Phase Shift Keying Modulation (DPSK) and Demodulation (DPSK)
        static Array<std::complex<double>> dpskmod(Array<double> data, int16_t M, double phase_offset);
        static Array<double> dpskdemod(Array<std::complex<double>> modulated_signal, int16_t M);

        // Pulse Amplitude Modulation (PAM) and Demodulation (PAM)
        static Array<std::complex<double>> pammod(Array<double> data, int16_t M, double phase_offset);
        static Array<double> pamdemod(Array<std::complex<double>> modulated_signal, int16_t M);

        // Pulse Width Modulation (PWM) and Demodulation (PWM)
        static Array<std::complex<double>> pwmmod(Array<double> data, int16_t M, double phase_offset);
        static Array<double> pwmdemod(Array<std::complex<double>> modulated_signal, int16_t M);

        // Pulse Position Modulation (PPM) and Demodulation (PPM)
        static Array<std::complex<double>> ppmmod(Array<double> data, int16_t M, double phase_offset);
        static Array<double> ppmdemod(Array<std::complex<double>> modulated_signal, int16_t M);
    };

    

   

} // namespace Com
