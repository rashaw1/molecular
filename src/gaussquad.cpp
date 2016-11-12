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
void GCQuadrature::initGrid(int points, GCTYPE _t) {
	t = _t;
	
	int offset, n;
	double nd;
	
	// Determine the nearest power of 2 <= points
	int power = (int) floor(log(points+1)/log(2));
	offset = pow(2, power);
	
	// Determine the integration parameters
	if (t == PS92) {
		order = offset - 1;
		N = order;
		M = (order - 1)/2;
		n = 1; 
		nd = n + 1.0;
		offset /= 2;
	} else if (t == PS93) {		
		order = 3*offset - 1;
		N = points+1;
		M = order/2;
		n = 3;
	}
	
	offset_fixed = offset; 
	
	// Allocate arrays and set parameters
	x = new double[order];
	w = new double[order];	
	start = 0;
	end = order - 1;
	
	// Determine weights and abscissae
	// Middle of interval
	x[M] = 0.0;
	w[M] = 1.0;
	// Common temporary variables
	double C1, S1, C0, S0, s, c, t;
	double opi2 = 2.0/M_PI;
	double opi23 = opi2 / 3.0;
	
	if (t == PS92) { // Following equations in PS92
		S0 = 1.0;
		C0 = 0.0;
		while (n <= M) {
			C1 = C0;
			S1 = S0;
			C0 = sqrt((1+C1)/2);
			S0 = S1/(2*C0);
			s = S0;
			c = C0;
		
			// Use recurrence relations
			offset /= 2;
			for (int i = 1; i <= n; i+=2) {
				t = 1 + opi23 * (3+2*s*s)*s*c - i/nd;
				int idx = i*offset - 1;
				x[order-idx-1] = t;
				x[idx] = -t;
				w[order - idx - 1] = w[idx] = s*s*s*s;
				t = s;
				s = s*C1 + c*S1;
				c = c*C1 - t*S1;
			}
			n = 2*n + 1;
			nd = n+1.0;
		}
	} else if (t == PS93) { // Following equations in PS93
		C0 = sin(M_PI/3);
		S0 = 0.5;
		C1 = S0;
		S1 = C0;
		c = cos(M_PI/3);
		s = C0;
		double s2 = s*s;
		
		t = (n - 2.0)/((double)n) + opi2 * (1 + 2.0*s2/3.0)*c*s;
		x[offset-1] = -t;
		x[order - offset] = t;
		w[order - offset] = w[offset-1] = s2*s2;
		
		while ((4*n/3 - 1) <= order) {
			c = C0;
			s = S0;
			offset /= 2;
			
			for (int i = 1; i < n; i += 2) {
				s2 = s*s;
				int idx = i*offset-1;
				t = 1 + opi23*s*c*(3 + 2*s2) - ((double)i)/((double)n);
				x[idx] = -t;
				x[order - idx - 1] = t;
				x[order - idx - 1] = w[idx] = s2*s2;
				
				t = s;
				s = s*C1 + c*S1;
				c = c*C1 - t*S1;
			}
			
			n *= 2;
			C1 = C0;
			S1 = S0;
			C0 = sqrt((1+C0)/2.0);
			S0 = S0/(2*C0);
		}
	}
	
	std::cout << order << " " << offset_fixed << " " << M << " " << start << " " << end << "\n";
	for (int q = 0; q < order; q++) std::cout << x[q] << " " << w[q] << "\n";
}

// Perform the GC integration on the function f
int GCQuadrature::integrate(std::function<double(double)> &f, const double tolerance) {
	int retval = -1; // 0 for converged, -1 for not converged
	
	// Initialise parameters
	int n = t == PS92 ? 1 : 3;
	double nd = n + 1.0;
	double e, T, q, p;
	int idx, i, cnt;
	int offset = offset_fixed;
	
	// PS92 Case
	if (t == PS92) {
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
				if (order - idx - 1 <= end) {
					T += w[order-idx-1]*f(x[order-idx-1]);
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
	} else if (t == PS93) {
		i, cnt, idx = offset - 1;
		p = w[M] * f(0.0);
		q = w[idx] * f(x[idx]) * f(x[order-offset]);
		I = p + q;
		offset /= 2;
		
		// Integrate
		int j = 0;
	    while((2*n*(1-j) + j*4*n/3 - 1) <= order) {
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
					if (order - idx - 1 <= end) {
						T += w[order-idx-1] * f(x[order-idx-1]);
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
// to the interval[0, infty), or [rmin, rmax]. We do this by the logarithmic transformation from KK
// or the linear mapping of FM06, respectively. 
void GCQuadrature::transformGrid(GCTYPE _t, double zeta_P, double P) {
	switch(_t) {
		case KK: {
			transformKK();
			break;
		}
		
		case FM06: {
			transformFM06(zeta_P, P);
			break;
		}
	}
}  

void GCQuadrature::transformKK() {
	double ln2 = log(2.0);
	double xt;
	
	for (int i = 0; i < order; i++) {
		xt = 1.0 - log(1.0-x[i])/ln2;
		w[i] = w[i]/(ln2 * (1.0 - x[i]));
		x[i] = xt;
	}
}

void GCQuadrature::transformFM06(double z, double P) {
	double sigma = 1.0 / sqrt(z);
	double t = P - 7.0 * sigma;
	double rmin = t > 0.0 ? t : 0.0;
	double rmax = P + 9.0*sigma;
	double i1 = 0.5*(rmax - rmin);
	double i2 = 0.5*(rmax + rmin);
	
	for (int i = 0; i < order; i++) {
		x[i] = i1*x[i] + i2;
		w[i] *= i1;
	}
}