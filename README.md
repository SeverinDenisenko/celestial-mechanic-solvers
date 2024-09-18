# adams-celestial-mechanic
 Adams solver for celestial mechanic applications

## What is it?

This is an implementation of two popular ODE solvers: the Runge-Kutta 4th order and Adams extrapolation methods (Linear multistep method). The main program solves a simple two-body problem as an example.

## Motivation

The Adams method of n'th order has n'th order of approximation, which makes it a flexible choice for celestial mechanics applications. The Adams method itself depends on other integrators to provide initial values for a multistep procedure. For this case was chosen Runge-Kutta 4th, but this can be changed.
