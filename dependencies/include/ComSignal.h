#pragma once
#include "Array.h"

namespace Com
{
    template <typename NUMERIC>
    class ComSignal
    {
    private:
        double Fs;
        Array<NUMERIC> samplePoints;
        
    public:
        ComSignal randi(double range, double n);
        ~ComSignal() = default;
    };
    
    

} // namespace Com
