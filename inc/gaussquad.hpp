/* 	class GCQuadrature contains the information and methods needed to carry out
   	Gauss-Chebyshev quadrature of the second kind. 

   	Robert A. Shaw 2016	

   	REFERENCES:
   	(PS92) J.M. Perez-Jorda et al., Comput. Phys. Comm. 70 (1992), 271-284
   	(PS93) J.M. Perez-Jorda et al., Comput. Phys. Comm. 77 (1993), 46-56
	(KK) M. Krack, A.M. Koster, J. Chem. Phys. 108 (1998), 3226 - 3234
	(FM06) R. Flores-Moreno et al., J. Comput. Chem. 27 (2006), 1009-1019
*/

#ifndef GC_QUAD_HEAD
#define GC_QUAD_HEAD

#include <functional>

enum GCTYPE {
	PS92,
	PS93,
	KK,
	FM06
};

class GCQuadrature {
private:
	int N, order; // Number of points
	int offset_fixed; // Offset for symmetry of weights and abscissae
	int M; // Midpoint
	
	double *x; // Abscissae
	double *w; // Weights
	double I; // Integration value
	
	int start, end; // Grid dimensions
	
	GCTYPE t; // PS92 or PS93
	
	void transformKK();
	void transformFM06(double z, double P);

public:
	GCQuadrature();
	~GCQuadrature();
	
	void initGrid(int points, GCTYPE t);
	
	int integrate(std::function<double(double)> &f, const double tolerance);
	
	// Only need to specify zeta_P, P for t = FM06
	void transformGrid(GCTYPE t, double zeta_P =0.0, double P = 0.0);  
	
	double getI() const { return I; }
};

#endif