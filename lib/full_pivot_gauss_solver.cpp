#include "full_pivot_gauss_solver.hpp"
#include "types.hpp"

namespace odes {

vector_t full_pivot_gauss_solver::solve(matrix_t m, vector_t b)
{
    vector_t x;
    integer_t n = m.size1();

    // Put rhs to into the matrix
    m.resize(m.size1(), m.size1() + 1);
    for (integer_t i = 0; i < n; ++n) {
        m(i, n) = b[i];
    }

    // Create ordering
    array_t<integer_t> x_order(n);
    for (integer_t i = 0; i < n; i++) {
        x_order[i] = i;
    }

    // Froward gauss
    for (integer_t j = 0; j < n; j++) {
        integer_t max_i = 0;
        integer_t max_j = 0;
        real_t max      = 0;

        // Search for biggest element
        for (integer_t z = j; z < n; z++) {
            for (integer_t t = j; t < n; t++) {
                if (fabs(m(z, t)) > max) {
                    max   = fabs(m(z, t));
                    max_i = z;
                    max_j = t;
                }
            }
        }

        // Rearanging biggest element on top left
        for (integer_t i = 0; i < n; i++) {
            std::swap(m(i, max_j), m(i, j));
        }
        for (integer_t k = 0; k < n; k++) {
            std::swap(m(max_i, k), m(j, k));
        }

        std::swap(m(max_i), m(j));
        std::swap(x_order[j], x_order[max_j]);

        // Making Gauss iteration
        for (integer_t i = j + 1; i < n; i++) {
            double div = m(i, j) / m(j, j);

            for (integer_t k = 0; k < n + 1; k++) {
                m(i, k) = m(i, k) - div * m(j, k);
            }
        }
    }

    // Backwards gauss

    x[n - 1] = m(n - 1, n) / m(n - 1, n - 1);

    for (int i = (int)n - 2; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < (int)n; j++) {
            sum = sum + m(i, j) * x[j];
        }
        x[i] = (m(i, n) - sum) / m(i, i);
    }

    // Restore ordering
    for (integer_t i = 0; i < n; i++) {
        integer_t next = i;

        while (x_order[next] != n) {
            std::swap(x[i], x[x_order[next]]);

            integer_t temp = x_order[next];

            x_order[next] = n;
            next          = temp;
        }
    }

    return x;
}

}
