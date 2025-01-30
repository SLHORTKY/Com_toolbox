#include <iostream>
#include "Array.h"
#include "SignalMath.h"
#include "matplotlibcpp.h"
#include "Noise.h"

namespace plt = matplotlibcpp;
using namespace Com;

double func(double x, std::unordered_map<std::string, double> &params)
{
    return x * params["x1"] * params["x2"];
}

int main()
{
    Array<double> x = Array<double>::arange(1, 10, 0.5);

    std::cout << x.apply(SignalMath::pow, 2).toString() << std::endl;

    Array<std::complex<double>> myArray = {
        std::complex<double>(1.0, 2.0), // complex number 1 (real=1.0, imag=2.0)
        std::complex<double>(3.0, 4.0), // complex number 2 (real=3.0, imag=4.0)
        std::complex<double>(5.0, 6.0), // complex number 3 (real=5.0, imag=6.0)
        std::complex<double>(7.0, 8.0)  // complex number 4 (real=7.0, imag=8.0)
    };

    std::cout<< myArray.apply(SignalMath::real).toString() << std::endl;

    std::cout << myArray.apply(SignalMath::rotate,45).toString() << std::endl;

}