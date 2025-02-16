#include <iostream>
#include "Array.h"
#include "SignalMath.h"
#include "ComSignal.h"
#include "matplotlibcpp.h"
#include "Filter.h"
#include "Modulator.h"
#include "Noise.h"

namespace plt = matplotlibcpp;
using namespace Com;

int main()
{
    
    Array<double> data = Array<double>::ones(100) * 0.0;
    data(1,99,20) = 5;

    ComSignal signal = ComSignal(data, Array<double>::linespace(0,100,100));

    plt::stem(signal.Domain(), signal.Signal());
    plt::show();

    ComSignal h = ComSignal(Array<double>{0.1,0.2,0.3,-0.1,-0.2,-0.3},Array<double>{0,1,2,3,4,5});

    ComSignal result = signal.convolution(h);

    plt::stem(result.Domain(), result.Signal());
    plt::show();


}