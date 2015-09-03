/*
 *
 *   PURPOSE: To declare a number of mathematical utility functions needed by 
 *            various other libraries in the molecular suite.
 * 
 *   LIST OF FUNCTIONS:
 *           fact(int i) : calculates i factorial
 *           fact2(int i) : calculates i double factorial
 *           boys(x, mmax, mmin) : calculates the boys function F_m(x)
 *                                 for m in the range mmax to mmin
 *                                 mmin = 0 by default. Returns values as
 *                                 a vector.
 *           binom(int m, int n) : calculates the binomial coefficient (n m)
 *
 *   DATE          AUTHOR            CHANGES
 *   =======================================================================
 *   27/08/15      Robert Shaw       Original code.
 *   28/08/15      Robert Shaw       Added boys function calculator.
 *   02/09/15      Robert Shaw       Added binomial coeff. calculator.
 */

#ifndef MATHUTILHEADERDEF
#define MATHUTILHEADERDEF

// Declare forward dependencies
class Vector;

// Functions to calculate the factorial and double factorial of an integer i
unsigned long int fact(int i);
unsigned long int fact2(int i);

// Boys function calculator for the evaluation of molecular integrals over
// gaussian type basis functions. Returns a vector of F_m(x) for m in the range
// mmax to mmin (it uses downwards recursion).
Vector boys(double x, int mmax, int mmin = 0, double PRECISION=1e-14);

// Calculates the binomial coefficient (n m)T
unsigned long int binom(int n, int m); 

#endif