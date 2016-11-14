/* 	class GCQuadrature contains the information and methods needed to carry out
   	Gauss-Chebyshev quadrature of the second kind. 

   	Robert A. Shaw 2016	

   	REFERENCES:
   	(Perez92) J.M. Perez-Jorda et al., Comput. Phys. Comm. 70 (1992), 271-284
   	(Perez93) J.M. Perez-Jorda et al., Comput. Phys. Comm. 77 (1993), 46-56
	(Krack98) M. Krack, A.M. Koster, J. Chem. Phys. 108 (1998), 3226 - 3234
	(Flores06) R. Flores-Moreno et al., J. Comput. Chem. 27 (2006), 1009-1019
*/

#ifndef GC_QUAD_HEAD
#define GC_QUAD_HEAD

#include <functional>

enum GCTYPE {
	ONEPOINT, // Described in Perez92
	TWOPOINT // Described in Perez93
};

class GCQuadrature {
private:
	int maxN; // Maximum number of points to use
	int M; // index of midpoint
	
	double *x; // Abscissae
	double *w; // Weights
	double I; // Integration value
	
	int start, end; // For prescreening
	
	GCTYPE t;
	
	double sumTerms(std::function<double(double)> &f, int limit, int shift, int skip);

public:
	GCQuadrature();
	~GCQuadrature();
	
	void initGrid(int points, GCTYPE t);
	
	// Returns true if quadrature converged, false otherwise. 
	bool integrate(std::function<double(double)> &f, const double tolerance);
	
	void transformZeroInf(); // Transformation from [-1, 1] to [0, infty) from Krack98
	void transformRMinMax(double z, double p);  // Transfromation from [-1, 1] to [rmin, rmax] from Flores06
	
	double getI() const { return I; }
};

#endif