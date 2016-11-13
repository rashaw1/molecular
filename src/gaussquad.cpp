/* 
	Implements gaussquad.hpp
	Robert A. Shaw 2016
 */

#include "gaussquad.hpp"
#include <cmath>
#include <iostream>

// Constructor
GCQuadrature::GCQuadrature() {
	// Currently does nothing
}

// Destructor
GCQuadrature::~GCQuadrature() {
	delete[] x;
	delete[] w;
}

// Initialise the quadrature grid
// As described in both Perez92 and Perez93
void GCQuadrature::initGrid(int points, GCTYPE _t) {
	t = _t;
	
	// Initialise parameters for grid
	int p;
	if (t == ONEPOINT) { // Perez92 one point method
		// We need the number of points to be of the form
		// 2^p - 1 for some power p. 
		p = (int) floor(log(points + 1)/log(2));
		maxN = pow(2, p) - 1;
	} else if (t == TWOPOINT) { // Perez93 two point method
		// Here we need it instead to be of the form
		// 3 * 2^p - 1 for some p.
		p = (int) floor(log((points + 2)/3.0)/log(2));
		maxN = 3*pow(2, p) - 1;
	}
	M = (maxN-1)/2; // Midpoint
	start = 0;
	end = maxN - 1;
	
	// initialise arrays
	x = new double[maxN];
	w = new double[maxN];
	
	// At the midpoint, M, x[M] = 0 and w[M] = 1
	x[M] = 0.0; w[M] = 1.0;
	// The rest of the abscissae and weights are then given by:
	// z_i = i*Pi / (maxN + 1), s_i = sin(z_i), c_i = cos(z_i)
	// x_i = 1 + 2/(3*pi) * [ (3 + 2*s_i^2) * c_i * s_i - 3z_i]
	// (3(maxN + 1)/16) * w_i = s_i^4
	// We then note that s_(i+1) = c_1 s_i + s_1 c_i
	// and c_(i+1) = c_1 c_i - s_1 s_i
	// with z_(i+1) = z_i + z_1
	// Clearly s_(maxN + 1 - i) = s_i
	// and c_(maxN + 1 - i) = -c_i
	// Therefore x_(maxN + 1 - i) = -x_i
	// and w_(maxN + 1 - i) = w_i
	// Meaning that we only have to calculate half the weights and abscissae
	double z1 = M_PI / ((double)(maxN + 1));
	double c1 = cos(z1); double s1 = sin(z1);
	double zi, si, ci, zi1, si1, ci1; //z_i, s_i, c_i, z_(i+1), s_(i+1), c_(i+1)
	zi1 = z1; si1 = s1; ci1 = c1;
	double o23pi = 2.0 / (3.0 * M_PI); // Convenient
	double s2; //si * si
	for (int n = 0; n < M; n++) {
		// First update zi, si, ci
		zi = zi1;
		si = si1;
		ci = ci1;
		s2 = si * si;
		
		// Now determine the w and x values
		w[maxN - 1 - n] = w[n] = s2 * s2;
		x[n] = 1 + o23pi * ( (3.0 + 2.0 * s2) * ci * si - 3.0*zi );
		x[maxN - 1 - n] = x[n];
		x[n] = -x[n];
		
		// Then update zi1, si1, ci1
		zi1 = zi + z1;
		si1 = c1 * si + s1 * ci;
		ci1 = c1 * ci - s1 * si;
	}
	
	/*std::cout << maxN << " " << M << " " << start << " " << end << "\n";
	for (int q = 0; q < maxN; q++) std::cout << x[q] << " " << w[q] << "\n";*/
}

// Perform the GC integration on the function f
int GCQuadrature::integrate(std::function<double(double)> &f, const double tolerance) {
	int retval = -1; // 0 for converged, -1 for not converged
	
	// Initialise parameters
	int n = t == ONEPOINT ? 1 : 3;
	double nd = n + 1.0;
	double e, T, q, p;
	int idx, i, cnt;
	int offset = (int) pow(2, (int) floor(log(maxN)/log(2)));
	
	// Perez92 Case
	if (t == ONEPOINT) {
		// Single point integration would use midpoint, M
		I = w[M]*f(x[M]);
		p = I;
		
		// Do the integration using an n-point Gaussian approximation 
		while(n<=M) {
			// Convergence checking variables
			q = 2*p; // q is 2I_(n-2)
			p = 2*I; // p is 2I_(n-1)
			offset /= 2; // Offset shrinks as n increases
			cnt = 0; // Count how many points were evaluated
		
			// Add extra points to integration for this n
			// Using relatinship I_(2n+1) = 0.5I_n + extra points
			for (i = 1; i <= n; i+=2) {
				idx = i*offset-1;
				T = 0.0;
				if (idx >= start) {
					T += w[idx]*f(x[idx]);
					cnt++;
				}
				if (maxN - idx - 1 <= end) {
					T += w[maxN-idx-1]*f(x[maxN-idx-1]);
					cnt++;
				}
				I += T;
			}
			n = 2*n + 1;
			nd = n+1.0;
		
			// Check convergence
			e = I - p;
			if (0==cnt) continue; // If no points were evaluated for this n, skip convergence check 
			if(16 * e * e <= 3 * nd * fabs(I - q) * tolerance) { 
				I = 16.0 * I / (3.0 * nd);
				retval = 0;
				break;
			}
				
		}
	} else if (t == TWOPOINT) {
		i, cnt, idx = offset - 1;
		p = w[M] * f(0.0);
		q = w[idx] * f(x[idx]) * f(x[maxN-offset]);
		I = p + q;
		offset /= 2;
		
		// Integrate
		int j = 0;
	    while((2*n*(1-j) + j*4*n/3 - 1) <= maxN) {
			j = 1 - j;
			if (j == 0) offset /= 2;
			cnt = 0;
			for (i=1; i<n; i+=2) {
				if (3 * ((i+2*j)/3) >= i+j) {
					idx = i*offset - 1;
					T = 0.0;
					if (idx >= start) {
						T += w[idx] * f(x[idx]);
						cnt++;
					}
					if (maxN - idx - 1 <= end) {
						T += w[maxN-idx-1] * f(x[maxN-idx-1]);
						cnt++;
					}
					I += T;
				}
			}
			i = n;
			n *= (1+j);
			p += (1-j)*(I-q);
			if (0 < cnt) e = 16*fabs((1-j)*(q-3*p/2) + j*(I-2*q)) / (3.0*n);
			q = (1-j)*q + j*I;
			if (0 == cnt) continue;
			
			if (e < tolerance) {
				I = 16.0 * q / (3.0 * n);
				retval = 0;
				break;
			}
		}	
	}
	
	return retval;
}

// The GC integrations above are over the interval [-1, 1] and thus need to be transformed
// to the interval[0, infty), or [rmin, rmax]. We do this by the logarithmic transformation from Krack98
// or the linear mapping of Flores06, respectively.  
void GCQuadrature::transformZeroInf() {
	double ln2 = log(2.0);
	double xt;
	
	for (int i = 0; i < maxN; i++) {
		xt = 1.0 - log(1.0-x[i])/ln2;
		w[i] = w[i]/(ln2 * (1.0 - x[i]));
		x[i] = xt;
	}
}

void GCQuadrature::transformRMinMax(double z, double p) {
	double osz = 1.0 / sqrt(z);
	
	// Determine interval
	double rmin = p - 7.0 * osz;
	rmin = rmin > 0 ? rmin : 0.0;
	double rmax = p + 9.0 * osz;
	
	// Find the relative and absolute midpoints 
	double rmid = 0.5*(rmax - rmin); // Midpoint of interval relative to rmin
	double amid = rmid + rmin; // Midpoint of interval
	
	// Transform weights and abscissae by linearly transforming
	// both are scaled by amid, and the abscissae are translated by amid
	for (int i = 0; i < maxN; i++) {
		x[i] = rmid * x[i] + amid;
		w[i] *= rmid;
	}
}