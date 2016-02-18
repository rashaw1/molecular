/*
 *    PURPOSE: To implement pbf.hpp, defining class PBF, a primitive
 *             cartesian gaussian basis function.
 *   
 *    DATE            AUTHOR              CHANGES
 *    =====================================================================
 *    27/08/15        Robert Shaw         Original code.
 *
 */

#include "pbf.hpp"
#include <cmath>
#include "mathutil.hpp"
#include <iostream>

PBF::PBF(double e, int l1, int l2, int l3) : exponent(e), lx(l1), ly(l2), lz(l3)
{
	normalise();  
}

// Copy constructor
PBF::PBF(const PBF& other)
{
  exponent = other.exponent;
  lx = other.lx; ly = other.ly; lz = other.lz;
  norm = other.norm;
}

// Routines

// Calculate the normalisation constant for a primitive cartesian gaussian
void PBF::normalise()
{
  // The formula can be found in Taketa, Huzinaga, and O-ohata, Journal of
  // the Physical Society of Japan, Vol. 21, No. 11, Nov 1966:
  // Gaussian-Expansion Methods for Molecular Integrals
  norm = std::pow(2, 2*(lx+ly+lz) + 1.5);
  norm = norm*std::pow(exponent, lx+ly+lz+1.5);
  // Calculate double factorials
  norm = norm / ( (double) (fact2(2*lx-1) * fact2(2*ly-1) * fact2(2*lz-1)) );
  norm = norm / std::pow(M_PI, 1.5);
  norm = std::sqrt(norm);
}

// Overloaded operators
PBF& PBF::operator=(const PBF& other)
{
  // Assign attributes
  exponent = other.exponent;
  norm = other.norm;
  lx = other.lx; ly = other.ly; lz = other.lz;
  return *this;
}
