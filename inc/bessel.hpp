/* 	class BesselFunction defines the weights and abscissae required for 
   	calculating the modified spherical Bessel function of the first kind,
   	needed for ECP integrals.

   	By Robert A. Shaw 2016 

	REFERENCES:
	R. Flores-Moreno et al., J. Comput. Chem. 27 (2005), 1009
	L.E. McMurchie, E. Davidson, J. Comp. Phys. 44 (1981), 289 
*/

#ifndef BESSEL_FUNCTION_HEAD
#define BESSEL_FUNCTION_HEAD

class Vector;

const double SMALL = 1.0E-7;
const int TAYLOR_CUT = 5;

class BesselFunction 
{
private:
	int lMax; // Maximum angular momentum
	int N; // Number of abscissae
	int order; // Order to which series is expanded
	
	double **K; // Bessel function values
	double *C; // Coefficients for derivatives of Bessel function
	
public:
	BesselFunction(int lMax, int N, int order, const double accuracy);
	~BesselFunction();
	
	int tabulate(const double accuracy);
	void calculate(const double z, Vector &values);
};

#endif
