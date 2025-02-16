#include "Filter.h"

using namespace Com;

Array<double> Com::Filter::generateSquareFilter(int length)
{
    return Array<double>::ones(length);
}

Array<double> Com::Filter::generateRaisedCosineFilter(int length, double rollOffFactor)
{
    return Array<double>();
}

Array<double> Com::Filter::generateRootRaisedCosineFilter(int length, double rollOffFactor)
{
    return Array<double>();
}

Array<double> Com::Filter::generateGaussianFilter(int length, double bandwidth)
{
    return Array<double>();
}

Array<double> Com::Filter::generateChebyshevFilter(int order, double ripple, double cutoffFrequency)
{
    return Array<double>();
}

Array<double> Com::Filter::generateMatchedFilter(const Array<double> &inputSignal)
{
    return Array<double>();
}

Array<double> Com::Filter::generateFIRFilter(const Array<double> &coefficients)
{
    return Array<double>();
}

Array<double> Com::Filter::generateIIRFilter(const Array<double> &aCoeffs, const Array<double> &bCoeffs)
{
    return Array<double>();
}

Array<double> Com::Filter::applyLMSFilter(const Array<double> &inputSignal, const Array<double> &desiredSignal, double stepSize)
{
    return Array<double>();
}

Array<double> Com::Filter::generateCombFilter(int delay, double gain)
{
    return Array<double>();
}
