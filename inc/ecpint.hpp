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
	
	FiveIndex U; // USP to RSH transformation coefficients
	FiveIndex W; // Type 1 angular integrals
	SevenIndex omega; // Type 2 angular integrals
	
	// Calculate terms for U and W
	double calcG(int l, int m, std::vector<double> &fac) const;
	double calcH1(int i, int j, int l, int m, std::vector<double> &fac) const;
	double calcH2(int i, int j, int k, int m, std::vector<double> &fac) const;
	ThreeIndex uklm(int lam, int mu, std::vector<double> &fac) const;
	ThreeIndex Pijk(int maxI) const; 
	
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
	
	bool isZero(int k, int l, int lam, int mu, double tolerance) const; 
	bool isZero(int k, int l, int m, int lam, int mu, double tolerance) const;
	bool isZero(int k, int l, int m, int lam, int mu, int rho, int sigma, double tolerance) const;
	
};

class RadialIntegral
{
};

class ECPIntegral
{
};


#endif