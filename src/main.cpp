#include <iostream>
#include "Array.h"
#include "SignalMath.h"
#include "matplotlibcpp.h"
#include "Filter.h"
#include "Modulator.h"
#include "Noise.h"

namespace plt = matplotlibcpp;

using namespace Com;

void func(){
    Noise noise;
    Array<double> data = Array<double>::ones(10); 
    auto upsampled_data = data.upSample(8);      

    // Generate a Root Raised Cosine (RRC) filter for zero ISI
    Array<double> rcos = Filter::generateRootRaisedCosineFilter(0.35, 8, 5);

    // Apply the RRC filter to shape the signal
    auto shaped_data = upsampled_data.applyFilter(rcos);

    // Get time-domain indices
    auto n_shaped = shaped_data.centeredDomain();
    auto n_rcos = rcos.centeredDomain();
    
    plt::stem(n_rcos, rcos);

    auto n_data = upsampled_data.centeredDomain();
    plt::figure();
    plt::stem(n_data, upsampled_data);

    plt::figure();
    plt::stem(n_shaped, shaped_data);

    auto down_sampled_shaped_data = shaped_data.downSample(8);

    plt::figure();
    plt::stem(down_sampled_shaped_data.centeredDomain(), down_sampled_shaped_data);
    plt::show();
}

int main()
{
    Noise noise;
    double M = 8;
    double N = 10;

    Array<double> data = noise.randi(0, M-1, N);
    Array<std::complex<double>> modulated = Modulator::ppmmod(data, M, 0.0);

    Array<double> t = Array<double>::linespace(0, modulated.size(),modulated.size());

    plt::stem(t , Array<double>::ones(modulated.size())* 0.1);
    plt::stem(t + modulated.apply(SignalMath::real) , Array<double>::ones(modulated.size()));
    plt::show();
}