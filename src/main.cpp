#include <iostream>
#include "Array.h"
#include "SignalMath.h"
#include "matplotlibcpp.h"
#include "Filter.h"
#include "Modulator.h"
#include "Noise.h"
#include "Helper.h"

namespace plt = matplotlibcpp;

using namespace Com;

void func()
{
    Noise n(42);

    int number_of_symbols = 1000;
    int M = 128;

    Array<double> data = n.randi(0, M - 1, number_of_symbols);

    Array<std::complex<double>> constellation = Modulator::qammod(data, M, false);

    plt::scatter(constellation.real(), constellation.imag());
    plt::show();
}

int main()
{

    auto Fs = 1000.0;                            
    auto T = 1.0/Fs;                  
    auto L = 2000;             
    auto t = Array<double>::arange(0,L-1,1)*T;

    auto S = t.apply([](double t){return SignalMath::sin<double>(t*2*M_PI*100);})*0.7 + t.apply([](double t){return SignalMath::sin<double>(t*2*M_PI*20);}) + 0.8;

    auto s_f = S.fft();

    plt::plot(Fs/L*Array<double>::arange(-L/2,L/2-1,1),s_f.fftshift().apply(SignalMath::abs<std::complex<double>>).real());
    plt::show();

    return 0;
}