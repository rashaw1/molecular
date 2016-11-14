/* class ECP contains the contracted expansion of primitive GaussianECPs that define a particular ECP
   class GaussianECP is simply a data structure for the primitive gaussian parameters.  

	Robert A. Shaw 2016

	REFERENCES:
	L.E. McMurchie, E.R. Davidson, J. Comput. Phys. 44 (1981), 289 - 301
	R. Flores-Moreno et al., J. Comp. Chem. 27 (2006), 1009 - 1019
 */

#ifndef ECP_HEAD
#define ECP_HEAD

#include <vector>

// Object describing a Gaussian of angular momentum l of the form
// d r^n e^{-ax^2}
struct GaussianECP {
	int n, l;
	double a, d;
	
	GaussianECP();
	GaussianECP(int n, int l, int a, int d);
	GaussianECP(const GaussianECP& other);
};

class ECP {
private:
	std::vector<GaussianECP> gaussians; // All the primitives in the ECP expansion
	int N, L; // # of Gaussians and maximum angular momentum
	
public:
	ECP();
	
	void addPrimitive(int n, int l, int a, int d, bool needSort = true);
	
	void sort(); // Sort primitives according to angular momentum
	
	// Evaluate U_l(r)
	double evaluate(double r, int l);
	
};

#endif
