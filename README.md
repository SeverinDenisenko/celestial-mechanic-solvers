# Celestial Mechanic Solvers
 Implementation of some ODE solvers for celestial mechanic.

## What is it?

This is an implementation of some popular ODE solvers.
* Runge-Kutta 4th order
* Adams extrapolation method
* Adams interpolation method
* Predictor-corrector method in 2 variants

The main program solves a simple two-body problem as an example.

Additional algorithms include:
* Quadrature folmuas
* * Trapezoidal rule
* * Simpson rule
* * Boole rule
* Full pivot Gauss method
* Jacobian matrix estimation
* * Second order derivatives
* * Fourth order derivatives
* Newton root finder

Some tests also can be used as examples.

## Details

All algorithms provide modern OO-design oriented interfaces.

This project uses the multiple-precision floating-point computation library MPFR (https://www.mpfr.org) and Boost uBLAS (https://github.com/boostorg/ublas). Everything is flexible enough to return to basic floating-point numbers.
