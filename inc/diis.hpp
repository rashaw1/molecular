/*
 *     PURPOSE: defines and implements a simple error class, containing
 *              an error code, and a message, and a way of printing
 *              to the primary output stream.
 * 
 *     DATE             AUTHOR                CHANGES
 *   =======================================================================
 *     14/08/15         Robert Shaw           Original code
 */

#ifndef DIISHEADERDEF
#define DIISHEADERDEF

// Includes
#include <vector> 
#include "mvector.hpp"
#include <Eigen/Dense>

class DIISEngine
{
private:
	int maxDiis;
	bool useDiis;
	std::vector<Vector> errs;
	double damping_factor;
	
	Eigen::MatrixXd lastB;
public:
	DIISEngine();
	void init(int _maxDiis, bool _useDiis, double _damp = 0.02);
	
	void use(bool on) { useDiis = on; }
	
	Vector compute(std::vector<Vector> &errors);
	Vector solve(Eigen::MatrixXd &B, int start);
};

#endif
