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

void func2 (){
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

int main()
{
    // Define radii and points per ring for 16-APSK
    const Array<double> radii = {1.0, 2.5};
    const Array<std::size_t> points_per_ring = {4, 12};

    // Generate LUT for 16-APSK constellation
    Array<std::complex<double>> apsk_lut = Com::generate_apsk_lut( points_per_ring, radii);

    // Print LUT points
    std::cout << "APSK LUT points:" << std::endl;
    for (std::size_t i = 0; i < apsk_lut.size(); ++i) {
        std::cout << i << ": " << apsk_lut[i] << std::endl;
    }

    // Sample data symbols to modulate (must be < 16)
    Array<double> data_symbols = {0, 5, 10, 15};

    // Map symbols to complex constellation points
    Array<std::complex<double>> modulated_signal = Modulator::apskmod(data_symbols, 16, radii, apsk_lut);

    // Print modulated symbols
    std::cout << "\nModulated signal:" << std::endl;
    for (auto& c : modulated_signal) {
        std::cout << c << std::endl;
    }

    return 0;
}