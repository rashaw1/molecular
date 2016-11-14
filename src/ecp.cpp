/* Implements ecp.hpp 

   Robert A. Shaw 2016
*/

#include "ecp.hpp"

#include <cmath>
#include <iostream>
#include <algorithm>

// GaussianECP constructor and copy constructor
GaussianECP::GaussianECP() : n(0), l(0), a(0), d(0) {}
GaussianECP::GaussianECP(int _n, int _l, int _a, int _d) : n(_n), l(_l), a(_a), d(_d) {}
GaussianECP::GaussianECP(const GaussianECP& other) : n(other.n), l(other.l), a(other.a), d(other.d) {}


// class ECP

ECP::ECP() : N(0), L(-1) {}

void ECP::addPrimitive(int n, int l, int a, int d, bool needSort) {
	GaussianECP newEcp(n, l, a, d);
	gaussians.push_back(newEcp);
	N++;
	L = l > L ? l : L;
	if (needSort) sort();
}

void ECP::sort() {
	std::sort(gaussians.begin(), gaussians.end(),
	 [&] (const GaussianECP& g1, const GaussianECP& g2) {return (g1.l < g2.l);});
}

// Evaluate U_l(r), assuming that gaussians sorted by angular momentum
double ECP::evaluate(double r, int l) {
	double value = 0.0;
	int am = 0;
	int i = 0;
	double r2 = r*r;
	while (am <= l && i < N) {
		am = gaussians[i].l;
		value += pow(r, gaussians[i].n) * gaussians[i].d * exp(-gaussians[i].a * r2);
		i++;
	} 
	return value; 
}

