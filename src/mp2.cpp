
#include "mp2.hpp"
#include "logger.hpp"
#include "mvector.hpp"
#include "matrix.hpp"
#include "molecule.hpp"
#include "integrals.hpp"
#include "error.hpp"
#include <iostream>
#include <thread>

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
	// Multithread
	int nthreads = focker.getMolecule().getLog().getNThreads();

	std::vector<Tensor4> moTemps(nthreads);
	std::vector<std::thread> thrds(nthreads);
	std::vector<int> startPoints(nthreads);
	std::vector<int> endPoints(nthreads);

	int spacing = N/nthreads;
	spacing = (double(N)/double(nthreads) - spacing > 0.5 ? spacing+1 : spacing);
    startPoints[0] = 0;
	endPoints[0] = spacing;
	
	for (int i = 0; i < nthreads; i++){
		moTemps[i] = Tensor4(endPoints[i]-startPoints[i], N, N, N, 0.0);
		thrds[i] = std::thread(&MP2::transformThread, *this, startPoints[i], endPoints[i], std::ref(moTemps[i]));
		if (i < nthreads - 1) { 
			startPoints[i+1] = endPoints[i];
			endPoints[i+1] = (i == (nthreads - 2) ? N : endPoints[i] + spacing);
		}
	}
	for (int i = 0; i < nthreads; i++) {
		thrds[i].join();

		// Copy in temporary matrix to moInts
		for (int p = startPoints[i]; p < endPoints[i]; p++){
			for (int q = 0; q < N; q++){
				for (int r = 0; r < N; r++){
					for (int s = 0; s < N; s++){
						moInts(p, q, r, s) = moTemps[i](p-startPoints[i], q, r, s);
					}
				}
			}
		}
	}
}

void MP2::transformThread(int start, int end, Tensor4& moTemp)
{
	Matrix& C = focker.getCP();
	IntegralEngine& aoInts = focker.getIntegrals();

	int offset = end - start;
	Tensor4 temp1(offset, N, N, N, 0.0);
	Tensor4 temp2(offset, N, N, N, 0.0);
	Tensor4 temp3(offset, N, N, N, 0.0);

	// Transform as four 'quarter-transforms'
	for (int p = start; p < end; p++){
		
		for (int mu = 0; mu < N; mu++){
			for (int a = 0; a < N; a++){
				for (int b = 0; b < N; b++){
					for (int c = 0; c < N; c++)
						temp1(p-start, a, b, c) += C(mu, p)*aoInts.getERI(mu, a, b, c);
				} // b
			} // a
		} // mu
		
		for (int q = 0; q < N; q++){
			for (int nu = 0; nu < N; nu++){
				for (int b = 0; b < N; b++){
					for (int c = 0; c < N; c++)
						temp2(p-start, q, b, c) += C(nu, q)*temp1(p-start, nu, b, c);
				} // b
			} // nu

			for (int r = 0; r < N; r++){
				for (int lam = 0; lam < N; lam++){
					for (int c = 0; c < N; c++)
						temp3(p-start, q, r, c) += C(lam, r)*temp2(p-start, q, lam, c);
				} // lam
				
				for (int s = 0; s < N; s++){
					for (int sig = 0; sig < N; sig++)
						moTemp(p-start, q, r, s) += C(sig, s)*temp3(p-start, q, r, sig);
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
					
