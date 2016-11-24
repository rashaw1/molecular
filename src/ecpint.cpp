/* Implements ecpint.hpp */

#include "matrix.hpp"
#include "ecpint.hpp"
#include "mathutil.hpp"
#include <iostream>

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

static int getIndex(int k, int l, int dim) { return k * dim + l; }
static int getIndex(int k, int l, int m, int dim1, int dim2) const { return k * dim1 * dim2 + l * dim2 + m; }
static int getIndex(int k, int l, int m, int n, int dim1, int dim2, int dim3) { return k * dim1 * dim2 * dim3 + l * dim2 * dim3 + m * dim3 + n; }

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


Matrix AngularIntegral::uklm(int lam, int mu, std::vector<double> &fac) const {
	Matrix values((lam+1)*(lam+1), 2); 
	 
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
		j = getIndex(k, l, lam+1);
		values(j, 0) = u;
		values(j, 1) = um;
	  }
	}
	return values;						
}


std::vector<double> AngularIntegral::pijk(int maxI) const {
	int dim = maxI+1;
	std::vector<double> values(dim * dim * dim);
	double pi4 = 4.0*M_PI;
	
	values[0] = pi4;
	for (int i = 1; i <= maxI; i++) {
		values[getIndex(i, 0, 0, dim, dim)] = pi4 / ((double) (2*i+1));
		
		for (int j = 1; j <= i; j++) {
			values[getIndex(i, j, 0, dim, dim)] = values[getIndex(i, j-1, 0, dim, dim)] * (2.0*j - 1.0) / (2.0 * ((double)(i + j)) + 1.0);
			
			for (int k = 1; k <= j; k++)
				values[getIndex(i, j, k, dim, dim)] = values[getIndex(i, j, k-1, dim, dim)] * (2.0*k - 1.0) / (2.0 * ((double)(i + j + k)) + 1.0);
			
		}
	}
	return values;
}

void AngularIntegral::makeU(std::vector<double> &fac) {
	int dim = maxLam + 1;
	
	int ix1, ix2;
	for (int lam = 0; lam <= maxLam; lam++) {
		for (int mu = 0; mu <= lam; mu++) {
			ix1 = getIndex(lam, mu, dim);
			Matrix Uij = uklm(lam, mu, fac);
			for (int i = 0; i <= lam; i++) {
				for (int j = 0; j <= lam; j++){
					ix2 = getIndex(i, j, dim);
					U(lam, mu, i, j, 0) = Uij(i, j, 0);
					U(lam, mu, i, j, 1) = Uij(i, j, 1);
				}
			}
		}
	}

	return U;

}



