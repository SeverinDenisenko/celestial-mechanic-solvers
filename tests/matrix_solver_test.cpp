#include <cstdlib>

#include "full_pivot_gauss_solver.hpp"
#include "types.hpp"

using namespace odes;

int main()
{
    matrix_t m(2, 2);
    vector_t a(2);
    vector_t b(2);
    vector_t c;
    m(0, 0) = 1;
    m(0, 1) = 0;
    m(1, 0) = 0;
    m(1, 1) = 1;
    a[0]    = 1;
    a[1]    = 1;
    b[0]    = 1;
    b[1]    = 1;

    odes::full_pivot_gauss_solver solver;
    c = solver.solve(m, a);

    if (norm(c - b) < 1e-10) {
        std::cout << "Good." << c << std::endl;
        exit(0);
    } else {
        std::cout << "Bad." << c << std::endl;
        exit(1);
    }
}
