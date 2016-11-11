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
	
	double *x; // Abscissae
	double *w; // Weights
	double I_value; // Integration value
	
	int start, end; // Grid dimensions

public:
	GCQuadrature();
	~GCQuadrature();
	
	void initGrid(int points, const double tolerance, GCTYPE t);
	
	int integrate(std::function<double(double, int, const void*)> &f, const void* params, GCTYPE t);
	
	// Only need to specify zeta_P, P for t = FM06
	void transformGrid(GCTYPE t, double zeta_P =0.0, double P = 0.0);  
};

#endif