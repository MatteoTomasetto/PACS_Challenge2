# Solvers to compute the zero of a function #

This program computes the zero of a function. In particular it considers f(t,y) = 0.5 - exp(pi * x).

This program allow the user to choose between different numerical methods in order to find the zero of a function:
- Bisection;
- RegulaFalsi;
- Brent;
- Secant;
- Newton;
- QuasiNewton.

The following parameters are taken in input from command line thanks to GetPot:
- method = name of the wanted method; (valid values: "Bisection", "RegulaFalsi", "Brent", "Secant", "Newton", "QuasiNewton");
- filename = name of the file with parameters written after the option -f or --file.

Example of execution: `./main method=Secant -f data`

By default: method = "Bisection" and filename = "data".

In this directory, `make` produces the executable which is just called `main`.

