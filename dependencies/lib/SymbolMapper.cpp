#include "SymbolMapper.h"

size_t Com::SymbolMapper::gray_encode(std::size_t n) {
    return n ^ (n >> 1);
}

size_t Com::SymbolMapper::gray_decode(std::size_t g) {
    std::size_t n = 0;
    for (; g; g >>= 1) {
        n ^= g;
    }
    return n;
}

std::vector<bool> Com::SymbolMapper::int_to_bits(std::size_t n, std::size_t bit_count) {
    if (bit_count == 0) return {};

    std::vector<bool> bits(bit_count);

    for (std::size_t i = 0; i < bit_count; ++i) {
        bits[bit_count - 1 - i] = (n >> i) & 1;
    }
    return bits;
}

std::size_t Com::SymbolMapper::bits_to_int(const std::vector<bool>& bits) {
    std::size_t value = 0;
    for (std::size_t i = 0; i < bits.size(); ++i) {
        value = (value << 1) | static_cast<std::size_t>(bits[i]);
    }
    return value;
}
