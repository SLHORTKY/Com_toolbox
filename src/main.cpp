#include <iostream>
#include "Array.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using namespace Com;

int main()
{
    Array x = Array::arange(1,20,1);

    std::cout<< x(1,10,1).toString() <<std::endl;
}