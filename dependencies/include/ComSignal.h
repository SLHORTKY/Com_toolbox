#pragma once
#include "Array.h"

namespace Com
{
    class ComSignal
    {
    private:
        double samplingFrequency;
        double samplingDuration;
        Array<std::complex<double>> samplePoints;
        
    public:
        
        ~ComSignal() = default;
    };
    
    

} // namespace Com
