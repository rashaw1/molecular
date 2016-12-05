#include "gshell.hpp"

GaussianShell::GaussianShell(double *A, int _l) : centerVec(A), l(_l) {}

void GaussianShell::addPrim(double e, double c) {
	exps.push_back(e);
	coeffs.push_back(c);
}