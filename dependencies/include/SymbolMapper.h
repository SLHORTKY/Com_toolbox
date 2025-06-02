#pragma once
#include <vector>
#include <cstddef>
#include <stdexcept>
#include <cmath>

namespace Com {
    class SymbolMapper {
    public:
        static std::size_t gray_encode(std::size_t n);
        static std::size_t gray_decode(std::size_t g);

        static std::vector<bool> int_to_bits(std::size_t n, std::size_t bit_count);
        static std::size_t bits_to_int(const std::vector<bool>& bits);
    };
}