#ifndef EXPONENTIALTIMEINTEGRATOR_H_
#define EXPONENTIALTIMEINTEGRATOR_H_

#include "boost/assign/std/vector.hpp"
using namespace boost::assign;

#include "TimeIntegrator.h"
#include "../utilities/BasicNames.h"
#include "../utilities/SemiParallelMatrix.h"

const int taylorSteps = 6;// Factor of the Taylor series; sets the number of iterations
const int arnoldiSize = 6; // Factor of the Arnoldi algorithm; sets the size of the orthonomal matrix V

namespace natrium {


/** @short Exponential time integration scheme for the solution of f' = L*f, as used in
 *         Uga etal. (2012) Spectral-element discontinuous Galerkin lattice Boltzmann simulation
 *         of flow past two cylinders in tandem with an exponential time integrator, CMWA 65 pp. 239-251
 */
template <class MATRIX, class VECTOR> class ExponentialTimeIntegrator : public TimeIntegrator <MATRIX, VECTOR> {

private:
	double factorial(int base);
	numeric_matrix makeIdentityMatrix();
	numeric_matrix makeMatrix(size_t m, size_t n);
	numeric_matrix taylorSeries (numeric_matrix base);
	numeric_matrix padeApproximation (numeric_matrix base);
	void phiFunction(const numeric_matrix &Hm, numeric_matrix &phiOne, numeric_matrix &phiTwo);

	numeric_matrix 	m_identityMatrix;
	numeric_matrix  m_Hm; 			// Hessenberg matrix (symmetric)
	numeric_matrix  m_H; 				// Hessenberg matrix
	numeric_matrix 	m_phiOne;
	numeric_matrix 	m_phiTwo;
	numeric_matrix 	m_phiExtended;
	VECTOR  	m_f;	 		//auxiliary vector for calculating f
	numeric_vector 	m_firstColumn; 	//vector for the first column of H_m (Symmetric Hessenberg matrix)

	VECTOR m_w; 	// Auxiliary vector for the Arnoldi algorithm
	VECTOR m_vj; // Auxiliary vector for the Arnoldi algorithm
	VECTOR m_vi; // Auxiliary vector for the Arnoldi algorithm

	SemiParallelMatrix<VECTOR> m_V;

	dealii::IndexSet getIndexSet (const MATRIX& m);

public:

	/// constructor
	ExponentialTimeIntegrator(double timeStepSize);


	ExponentialTimeIntegrator(double timeStepSize, size_t numberOfBlocks);

	/// destructor
	virtual ~ExponentialTimeIntegrator(){};

	virtual double step(VECTOR& f, const MATRIX& systemMatrix, const VECTOR& systemVector, double t = 0, double dt = 0);


};



} /* namespace natrium */
#endif /* EXPONENTIALTIMEINTEGRATOR_H_ */
