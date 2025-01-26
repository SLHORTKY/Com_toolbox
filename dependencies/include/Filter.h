#pragma once 
#include "ComSignal.h"

namespace Com
{
    class Filter : protected ComSignal
    {
    private:
        
    public:
        Filter(/* args */);
        ~Filter() = default;
    };
       
} // namespace Com
