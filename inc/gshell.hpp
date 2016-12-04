#ifndef GSHELL_HEAD
#define GSHELL_HEAD

#include <vector>

class GaussianShell {
private:
	std::vector<double> exps;
	std::vector<double> coeffs;
	double* centerVec;
	
public:
	GaussianShell(double* A);
	void addPrim(double exp, double c);
	
	int nprimitive() const { return exps.size(); }
	int ncartesian() const { return 6; }
	double* center() const { return centerVec; };
	double exp(int i) const { return exps[i]; }
	double coef(int i) const { return coeffs[i]; }
	int am() const { return 2; }
};

#endif