#pragma once
#include "ComSignal.h"

namespace Com
{
    enum class NOISE_TYPE
    {
    };

    class Noise : protected ComSignal
    {
    private:
        NOISE_TYPE type;

    public:
        Noise(/* args */);
        ~Noise() = default;
    };

} // namespace Com
