#pragma once

namespace Com
{
    enum class MODULATION_TYPE
    {
        AM,    // Amplitude Modulation: The amplitude of the carrier signal is varied in proportion to the message signal.
        FM,    // Frequency Modulation: The frequency of the carrier signal is varied according to the message signal.
        PM,    // Phase Modulation: The phase of the carrier signal is varied according to the message signal.
        QAM,   // Quadrature Amplitude Modulation: Combines both amplitude and phase modulation, allowing multiple bits per symbol.
        PAM,   // Pulse Amplitude Modulation: The amplitude of discrete pulses is varied according to the message signal.
        PWM,   // Pulse Width Modulation: The width of discrete pulses is varied in accordance with the message signal.
        PCM,   // Pulse Code Modulation: Analog signal is sampled and converted into digital code.
        DPSK,  // Differential Phase Shift Keying: A form of phase modulation where the phase difference between consecutive symbols represents the data.
        BPSK,  // Binary Phase Shift Keying: A phase modulation technique where binary 1 and 0 are represented by 180° phase shifts.
        QPSK,  // Quadrature Phase Shift Keying: A form of phase modulation where two bits are represented by four different phase states.
        FSK,   // Frequency Shift Keying: A modulation scheme where different frequencies represent binary 0 and 1.
        ASK,   // Amplitude Shift Keying: A modulation technique where binary data is represented by varying the amplitude of the carrier.
        OFDM,  // Orthogonal Frequency Division Multiplexing: A method of transmitting multiple signals using orthogonal subcarriers to maximize bandwidth efficiency.
        CPM,   // Continuous Phase Modulation: A phase modulation technique where the phase is continuously altered, improving spectral efficiency.
        BFSK,  // Binary Frequency Shift Keying: A special case of FSK using only two frequencies to represent binary 0 and 1.
        M_PSK, // M-ary Phase Shift Keying: A form of phase modulation using M different phase states to encode multiple bits per symbol.
        M_FSK, // M-ary Frequency Shift Keying: A form of FSK using M distinct frequencies to represent M different symbols.
        GMSK,  // Gaussian Minimum Shift Keying: A type of MSK with Gaussian filtering to reduce sideband power, used in GSM.
        FHSS,  // Frequency Hopping Spread Spectrum: A spread spectrum technique where the carrier hops between predefined frequencies to avoid interference.
        DSSS   // Direct Sequence Spread Spectrum: A spread spectrum technique where data is multiplied by a pseudorandom noise sequence to spread the signal over a wider bandwidth.
    };

    class Modulator
    {
    private:
        MODULATION_TYPE type;

    public:
        Modulator(/* args */);
        ~Modulator() = default;
    };

    Modulator::Modulator(/* args */)
    {
    }

} // namespace Com
