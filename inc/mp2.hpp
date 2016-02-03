#ifndef MP2HEADERDEF
#define MP2HEADERDEF

#include "tensor4.hpp"
#include "fock.hpp"

class MP2
{
private:
	int N, nocc;
	double energy;
	Tensor4 moInts;
	Fock& focker;
public:
	MP2(Fock& _focker);
	void transformIntegrals();
	void calculateEnergy();
	double getEnergy() const { return energy; }
};

#endif 
