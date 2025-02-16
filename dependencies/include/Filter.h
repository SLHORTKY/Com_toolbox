#pragma once 
#include "Array.h"

namespace Com
{
    class Filter 
    {
        
    public:
       // 1. Square Filter
    static Array<double> generateSquareFilter(int length);

    // 2. Raised Cosine Filter
    static Array<double> generateRaisedCosineFilter(int length, double rollOffFactor);

    // 3. Root Raised Cosine Filter
    static Array<double> generateRootRaisedCosineFilter(int length, double rollOffFactor);

    // 4. Gaussian Filter
    static Array<double> generateGaussianFilter(int length, double bandwidth);

    // 5. Butter worth Filter
    static std::array<double, 5> generateButterwortFilter(int order, double cutoffFrequency);

    // 6. Chebyshev Filter (Type I)
    static Array<double> generateChebyshevFilter(int order, double ripple, double cutoffFrequency);

    // 7. Matched Filter
    static Array<double> generateMatchedFilter(const Array<double>& inputSignal);

    // 8. FIR Filter
    static Array<double> generateFIRFilter(const Array<double>& coefficients);

    // 9. IIR Filter
    static Array<double> generateIIRFilter(const Array<double>& aCoeffs, const Array<double>& bCoeffs);

    // 10. Adaptive LMS Filter
    static Array<double> applyLMSFilter(const Array<double>& inputSignal, const Array<double>& desiredSignal, double stepSize);

    // 11. Comb Filter
    static Array<double> generateCombFilter(int delay, double gain);
    };
       
} // namespace Com
