/* 
	Implements gaussquad.hpp
	Robert A. Shaw 2016
 */

#include "gaussquad.hpp"
#include <cmath>

// Constructor
GCQuadrature::GCQuadrature() {
	// Currently does nothing
}

// Destructor
GCQuadrature::~GCQuadrature() {
	delete[] x;
	delete[] w;
}

// Initialise the quadrature grid
void GCQuadrature::initGrid(int points, const double tolerance, GCTYPE t) {
	// Determine the maximum order
	//order = pow(2, floor)
	
	// PS92 type
	
}