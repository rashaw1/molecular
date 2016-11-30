#include "gshell.hpp"

GaussianShell::GaussianShell(double *A) {
	centerVec = A;
}

void GaussianShell::addPrim(double e, double c) {
	exps.push_back(e);
	coeffs.push_back(c);
}