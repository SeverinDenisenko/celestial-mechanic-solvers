#include "utils.hpp"

namespace odes {

integer_t factorial(integer_t x)
{
    size_t res = 1;

    for (size_t i = 2; i <= x; ++i) {
        res *= i;
    }

    return res;
}

}
