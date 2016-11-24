/* ECP Integrals */

#ifndef ECPINT_HEAD
#define ECPINT_HEAD

#include <vector>

class Matrix;
class GCQuadrature;

// Calculate real spherical harmonics Slm(theta, phi) for all l, m up to lmax
static Matrix realSphericalHarmonics(int lmax, double theta, double phi, std::vector<double> &fac, std::vector<double> &dfac);  
static int getIndex(int k, int l, int dim);
static int getIndex(int k, int l, int m, int dim1, int dim2);
static int getIndex(int k, int l, int m, int n, int dim1, int dim2, int dim3);

class AngularIntegral 
{
private: 
	int LB, LE; // Maximum angular momentum of basis and ECP, respectively
	int wDim, ijkDim, maxL; // Limits for w-integral indices, and lambda
	
	Matrix U; // USP to RSH transformation coefficients
	Matrix W; // Type 1 angular integrals
	Matrix omega; // Type 2 angular integrals
	
	// Calculate terms for U and W
	double calcG(int l, int m, std::vector<double> &fac) const;
	double calcH1(int i, int j, int l, int m, std::vector<double> &fac) const;
	double calcH2(int i, int j, int k, int m, std::vector<double> &fac) const;
	Matrix uklm(int lam, int mu, std::vector<double> &fac) const;
	std::vector<double> pijk(int maxI) const; 
	
	// Build the angular integrals
	void makeU(std::vector<double> &fac);
	void makeW(std::vector<double> &fac);
	void makeOmega();
	
public:
	
	AngularIntegral(); // Default constructor
	AngularIntegral(int LB, int LE); 
	void init(int LB, int LE);
	void compute();
	void clear();
	
	double getU(int k, int l, int lam, int mu) const;
	double getIntegral(int k, int l, int m, int lam, int mu) const; 
	double getIntegral(int k, int l, int m, int lam, int mu, int rho, int sigma) const;
	
	bool isZero(int k, int l, int lam, int mu) const; 
	bool isZero(int k, int l, int m, int lam, int mu) const;
	bool isZero(int k, int l, int m, int lam, int mu, int rho, int sigma) const;
	
};

class RadialIntegral
{
};

class ECPIntegral
{
};


#endif