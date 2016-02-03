
#include "mp2.hpp"
#include "logger.hpp"
#include "mvector.hpp"
#include "matrix.hpp"
#include "molecule.hpp"
#include "integrals.hpp"
#include "error.hpp"
#include <iostream>

// Constructor
MP2::MP2(Fock& _focker) : focker(_focker)
{
	N = focker.getDens().nrows();
	nocc = focker.getMolecule().getNel()/2;
	energy = 0.0;
	moInts.assign(N, N, N, N, 0.0);
}

// Integral transformation
void MP2::transformIntegrals()
{
	Matrix& C = focker.getCP();
	IntegralEngine& aoInts = focker.getIntegrals();

	Tensor4 temp1(N, N, N, N, 0.0);
	Tensor4 temp2(N, N, N, N, 0.0);
	Tensor4 temp3(N, N, N, N, 0.0);
	
	// Transform as four 'quarter-transforms'
	for (int p = 0; p < N; p++){
		
		for (int mu = 0; mu < N; mu++){
			for (int a = 0; a < N; a++){
				for (int b = 0; b < N; b++){
					for (int c = 0; c < N; c++)
						temp1(p, a, b, c) += C(mu, p)*aoInts.getERI(mu, a, b, c);
				} // b
			} // a
		} // mu
		
		for (int q = 0; q < N; q++){
			for (int nu = 0; nu < N; nu++){
				for (int b = 0; b < N; b++){
					for (int c = 0; c < N; c++)
						temp2(p, q, b, c) += C(nu, q)*temp1(p, nu, b, c);
				} // b
			} // nu

			for (int r = 0; r < N; r++){
				for (int lam = 0; lam < N; lam++){
					for (int c = 0; c < N; c++)
						temp3(p, q, r, c) += C(lam, r)*temp2(p, q, lam, c);
				} // lam
				
				for (int s = 0; s < N; s++){
					for (int sig = 0; sig < N; sig++)
						moInts(p, q, r, s) += C(sig, s)*temp3(p, q, r, sig);
				} // s
			} // r
		} // q
	} // p
}

// Determine the MP2 energy
void MP2::calculateEnergy()
{
	Vector& eps = focker.getEps();

	double ediff, etemp;
	energy = 0.0;
	for (int i = 0; i < nocc; i++){
		for (int j = 0; j < nocc; j++){
			
			for (int a = nocc; a < N; a++){	
				for (int b = nocc; b < N; b++){
					ediff = eps[i] + eps[j] - eps[a] - eps[b];
					etemp = moInts(i, a, j, b)*(2.0*moInts(a, i, b, j) - moInts(a, j, b, i));
					
					energy += etemp/ediff;
				} // b
			} // a
		} // j
	} // i
}
					
