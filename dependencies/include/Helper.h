#pragma once
#include "Array.h"
#include <iostream>
#include <cmath>

namespace Com
{
    Com::Array<std::complex<double>> generate_apsk_lut(const Array<std::size_t> &counts, const Array<double> &radii)
    {
        if (counts.size() != radii.size())
            throw std::invalid_argument("counts and radii must match");

        Array<std::complex<double>> lut;
        for (std::size_t ring = 0; ring < counts.size(); ++ring)
        {
            std::size_t n = counts[ring];
            double r = radii[ring];
            for (std::size_t k = 0; k < n; ++k)
            {
                double theta = 2 * M_PI * k / n;
                lut.push_back(std::polar(r, theta));
            }
        }
        return lut;
    }
} // namespace Com
