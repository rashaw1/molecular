/* ECP Integrals */

#ifndef ECPINT_HEAD
#define ECPINT_HEAD

#include <vector>

class Matrix;
class GCQuadrature;

// Calculate single and double factorials iteratively
static std::vector<double> facArray(int l);
static std::vector<double> dfacArray(int l);

// Calculate real spherical harmonics Slm(theta, phi) for all l, m up to lmax
static Matrix realSphericalHarmonics(int lmax, double theta, double phi, std::vector<double> &fac, std::vector<double> &dfac);  

struct ThreeIndex {
	int dims[3];
	Matrix data;
	double& operator()(int i, int j, int k);
	double operator()(int i, int j, int k) const;
	ThreeIndex();
	ThreeIndex(int dim1, int dim2, int dim3);
	ThreeIndex(const ThreeIndex &other);
};

struct FiveIndex {
	int dims[5];
	Matrix data;
	double & operator()(int i, int j, int k, int l, int m);
	double operator()(int i, int j, int k, int l, int m) const;
	FiveIndex();
	FiveIndex(int dim1, int dim2, int dim3, int dim4, int dim5);
	FiveIndex(const FiveIndex &other);
};

struct SevenIndex {
	int dims[7];
	Matrix data;
	double & operator()(int i, int j, int k, int l, int m, int n, int p);
	double operator()(int i, int j, int k, int l, int m, int n, int p) const;
	SevenIndex();
	SevenIndex(int dim1, int dim2, int dim3, int dim4, int dim5, int dim6, int dim7);
	SevenIndex(const SevenIndex &other);
};

class AngularIntegral 
{
private: 
	int LB, LE; // Maximum angular momentum of basis and ECP, respectively
	int wDim, maxL; // Limits for w-integral indices, and lambda
	
	FiveIndex W; // Type 1 angular integrals
	SevenIndex omega; // Type 2 angular integrals
	
	// Calculate terms for U and W
	double calcG(int l, int m, std::vector<double> &fac) const;
	double calcH1(int i, int j, int l, int m, std::vector<double> &fac) const;
	double calcH2(int i, int j, int k, int m, std::vector<double> &fac) const;
	ThreeIndex uklm(int lam, int mu, std::vector<double> &fac) const;
	ThreeIndex Pijk(int maxI) const; 
	
	// Build the angular integrals
	FiveIndex makeU(std::vector<double> &fac);
	void makeW(std::vector<double> &fac, FiveIndex &U);
	void makeOmega(FiveIndex &U);
	
public:
	
	AngularIntegral(); // Default constructor
	AngularIntegral(int LB, int LE); 
	void init(int LB, int LE);
	void compute();
	void clear();
	
	double getIntegral(int k, int l, int m, int lam, int mu) const; 
	double getIntegral(int k, int l, int m, int lam, int mu, int rho, int sigma) const;

	bool isZero(int k, int l, int m, int lam, int mu, double tolerance) const;
	bool isZero(int k, int l, int m, int lam, int mu, int rho, int sigma, double tolerance) const;
	
};

class RadialIntegral
{
private:
  //	GCQuadrature bigGrid;
  //    GCQuadrature smallGrid;
	
	double tolerance;
	
	double T1Integrand(double r, double *p);
	double T2Integrand(double r, double *p);
	double Q2Integrand(double r, double *p);

public:
	RadialIntegral();
	void init(double tol = 1e-12, int small = 128, int large = 1024);
	
	std::vector<double> type1(double *coeffs, double *exponents);
	std::vector<double> type2(double *coeffs, double *exponents);	
};

class ECPIntegral
{
};


#endif
