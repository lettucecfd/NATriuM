#ifndef EXPONENTIALTIMEINTEGRATOR_H_
#define EXPONENTIALTIMEINTEGRATOR_H_

#include "boost/assign/std/vector.hpp"
using namespace boost::assign;

#include "TimeIntegrator.h"
#include "../utilities/BasicNames.h"

const int taylorSteps = 29;// Factor of the Taylor series; sets the number of iterations
const int arnoldiSize = 5; // Factor of the Arnoldi algorithm; sets the size of the orthonomal matrix V

namespace natrium {

/** @short Exponential time integration scheme for the solution of f' = L*f, as used in
 *         Uga etal. (2012) Spectral-element discontinuous Galerkin lattice Boltzmann simulation
 *         of flow past two cylinders in tandem with an exponential time integrator, CMWA 65 pp. 239-251
 */
template <class MATRIX, class VECTOR> class ExponentialTimeIntegrator : public TimeIntegrator <MATRIX, VECTOR> {

private:
	numeric_matrix makeIdentityMatrix()
	{
		numeric_matrix identityM(arnoldiSize);
		for (int i=0;i<arnoldiSize;i++)
				for(int j=0;j<arnoldiSize;j++)
				{
					{
						if(i==j)
							identityM(i,j)=1;
						else
							identityM(i,j)=0;
					}
				}

		return identityM;
	}

	numeric_matrix makeMatrix(size_t m, size_t n)
	{
		numeric_matrix newMatrix(m,n);
		return newMatrix;
	}
	// Calculates the matrix exponential of the given matrix;
	numeric_matrix taylorSeries(numeric_matrix basic)
	{
		numeric_matrix exponential(arnoldiSize);  // Builds the identity matrix
		numeric_matrix factor(arnoldiSize);
		exponential.copy_from(identityMatrix);

			factor = basic;
			exponential.add(factor,1);

			for(int j=2;j<taylorSteps;j++) // Taylor series as matrix exponential
			{
				factor.mmult(factor,basic);
				factor /=j;
				exponential.add (factor, 1);
			}
		return exponential;
	}

	numeric_matrix 	identityMatrix;
	numeric_matrix  H_m; 			// Hessenberg matrix (symmetric)
	numeric_matrix  H; 				// Hessenberg matrix
	numeric_vector 	m_f;	 		//auxiliary vector for calculating f
	numeric_vector 	firstColumn; 	//vector for the first column of H_m (Symmetric Hessenberg matrix)

	VECTOR w; 	// Auxiliary vector for the Arnoldi algorithm
	VECTOR v_j; // Auxiliary vector for the Arnoldi algorithm
	VECTOR v_i; // Auxiliary vector for the Arnoldi algorithm

public:

	/// constructor
	ExponentialTimeIntegrator(double timeStepSize, size_t problemSize);


	ExponentialTimeIntegrator(double timeStepSize, size_t problemSize, size_t numberOfBlocks);

	/// destructor
	virtual ~ExponentialTimeIntegrator(){};

	virtual void step(VECTOR& vector,
				const MATRIX& systemMatrix, const VECTOR& systemVector);


};





} /* namespace natrium */
#endif /* EXPONENTIALTIMEINTEGRATOR_H_ */
