/* Implements ecpint.hpp */

#include "matrix.hpp"
#include "ecpint.hpp"
#include "ecp.hpp"
#include "mathutil.hpp"
#include "gshell.hpp"
#include <iostream>
#include <functional>
#include <cmath>

// Compute single and double factorials iteratively
static std::vector<double> facArray(int l) {
	std::vector<double> values(l+1, 0.0);
	if (l > -1) {
		values[0] = 1.0;
		for (int i = 1; i < l + 1; i++) values[i] = values[i-1]*i;
	}
	return values; 
}

static std::vector<double> dfacArray(int l) {
	std::vector<double> values(l+1, 0.0);
	if (l > -1) {
		values[0] = 1.0;
		if (l > 0) {
			values[1] = 1.0;
			for (int i = 2; i <= l; i++) values[i] = values[i-2] * i;
		}
	}
	return values;
}

// Compute all the real spherical harmonics Slm(theta, phi) for l,m up to lmax
static Matrix realSphericalHarmonics(int lmax, double theta, double phi, std::vector<double> &fac, std::vector<double> &dfac){
	Matrix rshValues(lmax+1, 2*lmax+1);
	
	if (lmax > 0) {
		// First calculate the associated Legendre polynomials, Plm(cos theta), using the recursion relation
		// (l-m)Plm = x(2l - 1)P{l-1}m - (l+m-1)P{l-2}m
		// along with the zeroth order term
		// Pmm = (-1)^m (2m-1)!!(1-x^2)^{m/2}
		double x = cos(theta);
		double x2 = x * x;
		double Plm[lmax+1][lmax+1]; 
		// First get all Pmm terms
		Plm[0][0] = 1.0;
		double sox2 = sqrt(1.0 - x2);
		double ox2m = 1.0;
		for (int m = 1; m <= lmax; m++) {
			ox2m *= -sox2;
			Plm[m][m] = ox2m * dfac[2*m-1];
		}
		
		// Then increment l for each m
		Plm[1][0] = x;
		Plm[0][1] = 0.0;
		for (int l = 2; l <= lmax; l++) {
			ox2m = x * (2*l - 1);
			for (int m = 0; m < l; m++) {
				Plm[l][m] = ox2m * Plm[l-1][m] - (l + m - 1)*Plm[l-2][m];
				Plm[l][m] /= ((double) (l -m));
			}
			Plm[l-1][l] = 0.0;
		}
		
		// Now we compute the spherical harmonics via
		// Slm(theta, phi) = Clm * Plm(cos(theta)) * cos(m * phi), m > 0
		// Sl{-m}(theta, phi) = Clm * Plm(cos(theta)) * sin(m * phi)
		// Sl0(theta, phi) = sqrt(2) * Cl0 * Pl0(cos(theta))
		// where Clm^2 = (2l + 1)*(l - m)! / (8*pi * (l+m)!)
		double osq4pi = 1.0 / sqrt(4.0 * M_PI); 
		for (int l = 0; l <= lmax; l++) {
			rshValues(l, l) = osq4pi * sqrt(2.0 * l + 1.0) * Plm[l][0];
			for (int m = 1; m <= l; m++) {
				ox2m = (2.0 * l + 1.0) * fac[l-m] / fac[l+m];
				ox2m = osq4pi * sqrt(2.0 * ox2m) * Plm[l][m];
				rshValues(l, l+m) = ox2m * cos(m * phi);
				rshValues(l, l-m) = ox2m * sin(m * phi);
			}
		}
		
	} else {
		rshValues(0, 0) = 1.0 / sqrt(4.0 * M_PI);
	}
		
	return rshValues;
}

ThreeIndex::ThreeIndex() { dims[0] = 0; dims[1] = 0; dims[2] = 0; }
ThreeIndex::ThreeIndex(const ThreeIndex &other) { 
	data = other.data;
	for (int n = 0; n < 3; n++) dims[n] = other.dims[n]; 
}
ThreeIndex::ThreeIndex(int dim1, int dim2, int dim3) {
	dims[0] = dim1; dims[1] = dim2; dims[2] = dim3;
	data.resize(dim1, dim2*dim3);
}
double& ThreeIndex::operator()(int i, int j, int k) { return data(i, j*dims[2]+k); }
double ThreeIndex::operator()(int i, int j, int k) const { return data(i, j*dims[2]+k); }

FiveIndex::FiveIndex() { dims[0] = 0; dims[1] = 0; dims[2] = 0; dims[3] = 0; dims[4] = 0; }
FiveIndex::FiveIndex(const FiveIndex &other) { 
	data = other.data;
	for (int n = 0; n < 5; n++) dims[n] = other.dims[n]; 
}
FiveIndex::FiveIndex(int dim1, int dim2, int dim3, int dim4, int dim5) {
	dims[0] = dim1; dims[1] = dim2; dims[2] = dim3; dims[3] = dim4; dims[4] = dim5;
	data.resize(dim1*dim2, dim3*dim4*dim5);
}
double& FiveIndex::operator()(int i, int j, int k, int l, int m) { return data(i*dims[1] + j, k*dims[3]*dims[4] + l*dims[4] + m); }
double FiveIndex::operator()(int i, int j, int k, int l, int m) const { return data(i*dims[1] + j, k*dims[3]*dims[4] + l*dims[4] + m); }

SevenIndex::SevenIndex() { dims[0] = 0; dims[1] = 0; dims[2] = 0; dims[3] = 0; dims[4] = 0; dims[5] = 0; dims[6] = 0; }
SevenIndex::SevenIndex(const SevenIndex &other) { 
	data = other.data;
	for (int n = 0; n < 7; n++) dims[n] = other.dims[n]; 
}
SevenIndex::SevenIndex(int dim1, int dim2, int dim3, int dim4, int dim5, int dim6, int dim7) {
	dims[0] = dim1; dims[1] = dim2; dims[2] = dim3; dims[3] = dim4; dims[4] = dim5; dims[5] = dim6; dims[6]=dim7;
	data.resize(dim1*dim2*dim3, dim4*dim5*dim6*dim7);
}
double& SevenIndex::operator()(int i, int j, int k, int l, int m, int n, int p) { return data(i*dims[1]*dims[2] + j*dims[2] + k, l*dims[4]*dims[5]*dims[6] + m*dims[5]*dims[6] + n*dims[6] + p); }
double SevenIndex::operator()(int i, int j, int k, int l, int m, int n, int p) const { return data(i*dims[1]*dims[2] + j*dims[2] + k, l*dims[4]*dims[5]*dims[6] + m*dims[5]*dims[6] + n*dims[6] + p); }

double AngularIntegral::calcG(int l, int m, std::vector<double> &fac) const {
	double value = 0.0;
	double value1 = pow(2.0, l) * fac[l];
	value1 = 1.0 / value1; 
	double value2 = (2.0 * l + 1) * fac[l - m] / (2.0 * M_PI * fac[l + m]);
	value2 = sqrt(value2); 
	value = value1 * value2;
	return value;
} 

double AngularIntegral::calcH1(int i, int j, int l, int m, std::vector<double> &fac) const {
	double value = 0.0; 
	if (j > - 1){ 
		value = fac[l]/(fac[j]*fac[l - i]*fac[i-j]);
		value *= (1 - 2*(i%2)) * fac[2*(l - i)] / (fac[l - m - 2*i]);
	}
	return value;
}

double AngularIntegral::calcH2(int i, int j, int k, int m, std::vector<double> &fac) const {
	double value = 0.0; 
	int ki2 = k - 2*i;
	if ( m >= ki2 && ki2 >= 0 ) {
		value = fac[j]*fac[m]/(fac[i] * fac[j-i] * fac[ki2] * fac[m-ki2]);
		int p = (m - k + 2*i)/2;
		value *= (1.0 - 2.0*(p%2));
	}
	return value;
}


ThreeIndex AngularIntegral::uklm(int lam, int mu, std::vector<double> &fac) const {
	ThreeIndex values(lam+1, lam+1, 2);
	 
  	double or2 = 1.0/sqrt(2.0);
  	double u = 0.0;
	double um = 0.0;
	double g = calcG(lam, mu, fac);

  	double u1, h1, h2;
  	int j;
  	for (int k = 0; k <= lam; k++) {
  	  for (int l = 0; l <= lam - k; l++) {
		u = um = 0.0;
	  	j = k + l - mu;
		if (j % 2 == 0) { 
			u1 = 0.0;
			j/=2;
			for (int i = j; i <= (lam - mu)/2; i++) u1 += calcH1(i, j, lam, mu, fac);
			
			u = g * u1;
			u1 = 0;
			for (int i = 0; i <= j; i++) u1 += calcH2(i, j, k, mu, fac);
			u *= u1;
			um = u;
			
			j = l % 2;
			u *= (1 - j);
			um *= j;
			if (mu == 0) {
				u *= or2;
				um = u;
			} 
		}
		values(k, l, 0) = u;
		values(k, l, 1) = um;
	  }
	}
	return values;						
}


ThreeIndex AngularIntegral::Pijk(int maxI) const {
	int dim = maxI+1;
	ThreeIndex values(dim, dim, dim);
	double pi4 = 4.0*M_PI;
	
	values(0, 0, 0) = pi4;
	for (int i = 1; i <= maxI; i++) {
		values(i, 0, 0) = pi4 / ((double) (2*i+1));
		
		for (int j = 1; j <= i; j++) {
			values(i, j, 0) = values(i, j-1, 0) * (2.0*j - 1.0) / (2.0 * ((double)(i + j)) + 1.0);
			
			for (int k = 1; k <= j; k++)
				values(i, j, k) = values(i, j, k-1) * (2.0*k - 1.0) / (2.0 * ((double)(i + j + k)) + 1.0);
			
		}
	}
	return values;
}

FiveIndex AngularIntegral::makeU(std::vector<double> &fac) {
	int dim = maxL + 1;

	FiveIndex values(dim, dim, dim, dim, 2);
	for (int lam = 0; lam <= maxL; lam++) {
		for (int mu = 0; mu <= lam; mu++) {
			ThreeIndex Uij = uklm(lam, mu, fac);
			for (int i = 0; i <= lam; i++) {
				for (int j = 0; j <= lam; j++){
					values(lam, mu, i, j, 0) = Uij(i, j, 0);
					values(lam, mu, i, j, 1) = Uij(i, j, 1);
				}
			}
		}
	}
	
	return values;
}

void AngularIntegral::makeW(std::vector<double> &fac, FiveIndex &U) {
	int LB2 = 2*LB;
	int dim = wDim;
	int maxI = dim/2 + LB;
	int maxLam = maxL;
	
	FiveIndex values{dim+1, dim+1, dim+1, maxLam+1, 2*(maxLam + 1)};
	ThreeIndex pijk = Pijk(maxI);
	
	int plam, pmu;
	double smu, w;
	std::vector<int> ix(3);
	for (int k = 0; k <= dim; k++) {	
		for (int l = 0; l <= dim; l++) {	
			for(int m = 0; m <= dim; m++) {
				plam = (k + l + m)%2;
				
				int limit = maxLam > k+l+m ? k+l+m : maxLam;
				for(int lam = plam; lam <= limit; lam += 2){
					smu = 1 - 2*(l%2);
					pmu = (k+l) % 2;
					
					for (int mu = pmu; mu <= lam; mu+=2) {
						w = 0.0;
						
						for (int i = 0; i <= lam; i++) {
							for (int j = 0; j <= lam - i; j++) {
								ix[0] = k+i;
								ix[1] = l+j;
								ix[2] = m + lam - i - j; 
								
								if (ix[0]%2 + ix[1]%2 + ix[2]%2 == 0){
									std::sort(ix.begin(), ix.end()); 
									w += U(lam, mu, i, j, (1 - (int)(smu))/2)*pijk(ix[2]/2, ix[1]/2, ix[0]/2);
								}
							}
						}
						
						values(k, l, m, lam, lam+(int)(smu*mu)) = w;
					}
				}	
			}	
		}	
	}
	W = values;
}

void AngularIntegral::makeOmega(FiveIndex &U) {
	
	int lamDim = LE + LB; 
	int muDim = 2*lamDim + 1;
	SevenIndex values{LB+1, LB+1, LB+1, lamDim+1, muDim+1, lamDim+1, muDim+1};
	
	double om_plus=0.0, om_minus=0.0;
	double wval; 
	for (int k = 0; k <= LB; k++) {
		for (int l = 0; l <= LB; l++) {
			for (int m = 0; m <= LB; m++) {
					
				for (int rho = 0; rho <= lamDim; rho++ ) {
					for (int sigma = -rho; sigma <= rho; sigma++) {
						
						for (int lam = 0; lam <= rho; lam++) {
							for (int mu = 0; mu <= lam; mu++) {
								
								om_plus = om_minus = 0.0;
								for (int i = 0; i<= lam; i++ ) {
									for (int j = 0; j <= lam - i; j++) {
										
										wval = W(k+i, l+j, m+lam-i-j, rho, rho+sigma);
										om_plus += U(lam, mu, i, j, 0) * wval;
										om_minus += U(lam, mu, i, j, 1) * wval;
										
									}
								}
								if (mu == 0) om_minus = om_plus;
								values(k, l, m, rho, sigma+rho, lam, lam+mu) = om_plus;
								values(k, l, m, lam, lam+mu, rho, sigma+rho) = om_plus;
								values(k, l, m, rho, sigma+rho, lam, lam-mu) = om_minus;
								values(k, l, m, lam, lam-mu, rho, sigma+rho) = om_minus;
								
							}
						}
						
					}
				}
					
			}
		}
	}
	
	omega = values;
}

AngularIntegral::AngularIntegral() { init(0, 0); }
AngularIntegral::AngularIntegral(int _LB, int _LE) { init(_LB, _LE); }
void AngularIntegral::init(int _LB, int _LE ) {
	LB = _LB;
	LE = _LE;
	wDim = 4*LB > 3*LB + LE ? 4*LB : 3*LB + LE;
	maxL = 2*LB > LB + LE ? 2*LB : LB+LE;
	
}

void AngularIntegral::compute() {
	std::vector<double> fac = facArray(wDim);
	
	FiveIndex U = makeU(fac);
	makeW(fac, U);
	makeOmega(U);
}

void AngularIntegral::clear() {}

double AngularIntegral::getIntegral(int k, int l, int m, int lam, int mu) const { return W(k, l, m, lam, lam+mu); }
double AngularIntegral::getIntegral(int k, int l, int m, int lam, int mu, int rho, int sigma) const { return omega(k, l, m, lam, lam+mu, rho, rho+sigma); }

bool AngularIntegral::isZero(int k, int l, int m, int lam, int mu, double tolerance) const {
	if (wDim > 0) return fabs(W(k, l, m, lam, lam+mu)) < tolerance;
	else return true;
}
bool AngularIntegral::isZero(int k, int l, int m, int lam, int mu, int rho, int sigma, double tolerance) const {
	if (wDim > 0) return fabs(omega(k, l, m, lam, lam+mu, rho, rho+sigma)) < tolerance;
	else return true;
}

//****************************************** RADIAL INTEGRAL *********************************************

RadialIntegral::RadialIntegral() {}

void RadialIntegral::init(int maxL, double tol, int small, int large) {
	bigGrid.initGrid(large, ONEPOINT);
	smallGrid.initGrid(small, TWOPOINT);
	smallGrid.transformZeroInf();
	
	bessie.init(maxL, 1600, 200, tol);
	
	tolerance = tol;
}

void RadialIntegral::buildBessel(double *r, int nr, int maxL, Matrix &values, double weight) {
	std::vector<double> besselValues;
	for (int i = 0; i < nr; i++) {
		bessie.calculate(weight * r[i], maxL, besselValues);
		for (int l = 0; l <= maxL; l++) values(l, i) = besselValues[l];
	}
}

double RadialIntegral::calcKij(double Na, double Nb, double zeta_a, double zeta_b, double *A, double *B) const {
	double muij = zeta_a * zeta_b / (zeta_a + zeta_b);
	double R[3] = {A[0] - B[0], A[1] - B[1], A[2] - B[2]};
	double R2 = R[0] * R[0] + R[1] * R[1] + R[2] * R[2];
	return Na * Nb * exp(muij * R2);
}

// Assumes that p is the pretabulated integrand at the abscissae
double RadialIntegral::integrand(double r, double *p, int ix) {
	return p[ix];
}

void RadialIntegral::buildParameters(GaussianShell &shellA, GaussianShell &shellB, double *A, double *B) {
	int npA = shellA.nprimitive();
	int npB = shellB.nprimitive();

	p.assign(npA, npB, 0.0);
	P.assign(npA, npB, 0.0);
	P2.assign(npA, npB, 0.0);
	K.assign(npA, npB, 0.0);

	double Pvec[3];
	double zetaA, zetaB;
	for (int a = 0; a < npA; a++) {
		zetaA = shellA.exp(a);
		
		for (int b = 0; b < npB; b++) {
			zetaB = shellB.exp(a);
			
			p(a, b) = zetaA + zetaB;
			for (int n = 0; n < 3; n++) {
				Pvec[n] = (zetaA * A[n] + zetaB * B[n])/p(a, b);
				P2(a, b) = Pvec[0]*Pvec[0] + Pvec[1]*Pvec[1] + Pvec[2]*Pvec[2];
				P(a, b) = sqrt(P2(a, b));
			}
			K(a, b) = calcKij(1.0, 1.0, zetaA, zetaB, A, B);
			
		}
	}
}

void RadialIntegral::buildU(ECP &U, int l, GaussianShell &shellA, GaussianShell &shellB, GCQuadrature &grid, double *Utab) {
	int gridSize = smallGrid.getN();
	double* gridPoints = smallGrid.getX();
	
	// Tabulate weighted ECP values
	int N = shellA.am() + shellB.am();
	double r;
	bool foundStart = false;
	for (int i = 0; i < gridSize; i++) {
		r = gridPoints[i];
		Utab[i] = pow(r, N+2) * U.evaluate(r, l);
		if(Utab[i] > tolerance && !foundStart) {
			 grid.start = i;
			 foundStart = true;
		}
		if(Utab[i] < tolerance && foundStart) {
			grid.end = i-1;
			foundStart = false;
		}
	}
}

int RadialIntegral::integrate(int maxL, int gridSize, Matrix &intValues, GCQuadrature &grid, std::vector<double> &values) {
	std::function<double(double, double*, int)> intgd = integrand; 
	values.assign(maxL+1, 0.0);
	int test;
	double params[gridSize];
	for (int i = 0; i < grid.start; i++) params[i] = 0.0;
	for (int i = grid.end+1; i < gridSize; i++) params[i] = 0.0;
	for (int l = 0; l <= maxL; l++) {
		for (int i = grid.start; i <= grid.end; i++) params[i] = intValues(l, i);
		test = grid.integrate(intgd, params, tolerance);
		values[l] = smallGrid.getI();
		if (test == 0) break;
	}
	return test;
}

void RadialIntegral::type1(int maxL, ECP &U, GaussianShell &shellA, GaussianShell &shellB, double *A, double *B, std::vector<double> &values) {
	int npA = shellA.nprimitive();
	int npB = shellB.nprimitive();
	
	buildParameters(shellA, shellB, A, B);
	
	int gridSize = smallGrid.getN();
	double* gridPoints = smallGrid.getX();
	
	// Reset grid starting points
	smallGrid.start = 0;
	smallGrid.end = gridSize;
	
	double Utab[gridSize];
	buildU(U, U.getL(), shellA, shellB, smallGrid, Utab);

	// Now pretabulate integrand
	Matrix intValues(maxL+1, gridSize, 0.0);
	// and bessel function
	Matrix besselValues(maxL+1, gridSize);
	// Calculate type1 integrals
	double da, db, val;
	// Tabulate integrand
	for (int a = 0; a < npA; a++) {
		da = shellA.coef(a);
		for (int b = 0; b < npB; b++) {
			db = shellB.coef(b);
			
			buildBessel(gridPoints, gridSize, maxL, besselValues, 2.0*p(a,b)*P(a,b));
			
			for (int i = smallGrid.start; i <= smallGrid.end; i++) {
				val = -p(a, b) * (gridPoints[i]*(gridPoints[i] - 2*P(a, b)) + P2(a, b));
				val = da * db * K(a, b) * exp(val);
			
				for (int l = 0; l <= maxL; l++)
					intValues(l, i) += val * besselValues(l, i);
				
				for (int l = 0; l <= maxL; l++)
					intValues(l, i) *= Utab[i];
			}
		}
	}
	
	// Calculate integrals
	int test = integrate(maxL, gridSize, intValues, smallGrid, values);
	if (test == 0) std::cout << "Failed to converge\n";
}

// F_a(lam, r) = sum_{i in a} d_i K_{lam}(2 zeta_a A r)*exp(-zeta_a(r - A)^2)
void RadialIntegral::buildF(GaussianShell &shell, double *Avec, int maxL, double *r, int nr, int start, int end, Matrix &F) {
	int np = shell.nprimitive();
	double A = Avec[0]*Avec[0] + Avec[1]*Avec[1] + Avec[2]*Avec[2];
	A = sqrt(A); 
		
	double weight, zeta, c;
	Matrix besselValues(maxL+1, nr, 0.0);
	
	F.assign(maxL+1, nr, 0.0);
	for (int a = 0; a < np; a++) {
		zeta = shell.exp(a);
		c = shell.coef(a);
		weight = 2.0 * zeta * A;
		
		buildBessel(r, nr, maxL, besselValues, weight);
		
		for (int i = start; i <= end; i++) {
			weight = r[i] - A;
			weight = c * exp(-zeta * weight * weight);
			
			for (int l = 0; l <= maxL; l++) 
				F(l, i) += weight * besselValues(l, i); 
		}
	}
}

void RadialIntegral::type2(int l, int maxL1, int maxL2, ECP &U, GaussianShell &shellA, GaussianShell &shellB, double *Avec, double *Bvec, Matrix &values) {
	int npA = shellA.nprimitive();
	int npB = shellB.nprimitive();
	
	buildParameters(shellA, shellB, Avec, Bvec);
	
	std::cout << "Built params\n";
	
	// Start with the small grid
	// Pretabulate U
	int gridSize = smallGrid.getN();
	double* gridPoints = smallGrid.getX();
	
	// Reset grid starting points
	smallGrid.start = 0;
	smallGrid.end = gridSize;
	
	double Utab[gridSize];
	buildU(U, l, shellA, shellB, smallGrid, Utab);
	
	// Build the F matrices
	Matrix Fa;
	Matrix Fb;
	buildF(shellA, Avec, maxL1, gridPoints, gridSize, smallGrid.start, smallGrid.end, Fa);
	buildF(shellB, Bvec, maxL2, gridPoints, gridSize, smallGrid.start, smallGrid.end, Fb);
	
	// Build the integrals
	Matrix intValues(maxL2+1, gridSize, 0.0);
	std::vector<int> tests(maxL1+1);
	std::vector<double> tempValues;
	bool failed = false;
	values.assign(maxL1+1, maxL2+1, 0.0);
	for (int l1 = 0; l1 <= maxL1; l1++) {
		for (int i = smallGrid.start; i <= smallGrid.end; i++) {
			for (int l2 = 0; l2 <= maxL2; l2++) 
				intValues(l2, i) = Utab[i] * Fa(l1, i) * Fb(l2, i);
		}
		tests[l1] = integrate(maxL2, gridSize, intValues, smallGrid, tempValues);
		failed = failed || (tests[l1] == 0);
		for (int l2 = 0; l2 <= maxL2; l2++) values(l1, l2) = tempValues[l2];
	}
	
	if (failed) {
		// Not converged, switch to big grid
		GCQuadrature newGrid;
		double zeta_a, zeta_b, c_a, c_b, weight, XA, XB;
		double A = Avec[0]*Avec[0] + Avec[1]*Avec[1] + Avec[2]*Avec[2];
		double B = Bvec[0]*Bvec[0] + Bvec[1]*Bvec[1] + Bvec[2]*Bvec[2];
		A = sqrt(A); B = sqrt(B);
		
		for (int l1 = 0; l1 <= maxL1; l1++) {
			if (tests[l1] == 0) { 
				for (int l2 = 0; l2 <= maxL2; l2++) values(l1, l2) = 0.0;
			
				for (int a = 0; a < npA; a++) {
					zeta_a = shellA.exp(a);
					c_a = shellA.coef(a);
				
					// Build bessel function values
					weight = 2.0 * zeta_a * A;
					buildBessel(gridPoints, gridSize, maxL2, Fa, weight);
					for (int i = 0; i < gridSize; i++) {
						XA = gridPoints[i] - A;
						XA = exp(-zeta_a * XA * XA);
						for (int l2 = 0; l2 <= maxL2; l2++) Fa(l2, i) *= XA;
					}
				
					// calculate exponential
					
					for (int b = 0; b < npB; b++) {
						zeta_b = shellB.exp(b);
						c_b = shellB.coef(b); 
					
						// Build bessel function values
						weight = 2.0 * zeta_b * B;
					
						// Set up grid
						newGrid = bigGrid;
						gridSize = newGrid.getN();
						gridPoints = newGrid.getX();
						newGrid.start = 0;
						newGrid.end = gridSize;
						newGrid.transformRMinMax(p(a,b), (zeta_a * A + zeta_b * B)/p(a, b));
				
						// Build the U tab
						buildU(U, l, shellA, shellB, newGrid, Utab);				
						intValues.assign(maxL2+1, gridSize, 0.0);
					
						// Build U and bessel
					
						buildBessel(gridPoints, gridSize, maxL2, Fb, weight); 
						for (int i = 0; i < gridSize; i++) {
							XB = gridPoints[i] - B;
							XB = exp(-zeta_b * XB * XB);
							for (int l2 = 0; l2 <= maxL2; l2++) {
								Fb(l2, i) *= XB;
								intValues(l2, i) = Utab[i] * Fa(l2, i) * Fb(l2, i);
							}		
						}
					
						integrate(maxL2, gridSize, intValues, newGrid, tempValues);
						for (int l2 = 0; l2 <= maxL2; l2++) values(l1, l2) += c_a*c_b*tempValues[l2];
					}
				}
			}
		}
	}
	
}


