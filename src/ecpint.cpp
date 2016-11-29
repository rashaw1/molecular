/* Implements ecpint.hpp */

#include "matrix.hpp"
#include "ecpint.hpp"
#include "mathutil.hpp"
#include <iostream>
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

