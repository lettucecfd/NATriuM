#ifndef EXPONENTIALTIMEINTEGRATOR_H_
#define EXPONENTIALTIMEINTEGRATOR_H_

#include "boost/assign/std/vector.hpp"
using namespace boost::assign;

#include "TimeIntegrator.h"
#include "../utilities/BasicNames.h"

const int taylorSteps = 6;// Factor of the Taylor series; sets the number of iterations
const int arnoldiSize = 6; // Factor of the Arnoldi algorithm; sets the size of the orthonomal matrix V

namespace natrium {

/** @short Exponential time integration scheme for the solution of f' = L*f, as used in
 *         Uga etal. (2012) Spectral-element discontinuous Galerkin lattice Boltzmann simulation
 *         of flow past two cylinders in tandem with an exponential time integrator, CMWA 65 pp. 239-251
 */
template <class MATRIX, class VECTOR> class ExponentialTimeIntegrator : public TimeIntegrator <MATRIX, VECTOR> {

private:
	double factorial(int base)
	{
		unsigned int factorial =1;
		for (int i=1;i<=base;i++)
		{
			factorial = factorial * i;
		}
		return factorial;
	}

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
							identityM(i,j)=0 ;
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
	numeric_matrix taylorSeries(numeric_matrix base)
	{
		numeric_matrix exponential(arnoldiSize);  // Builds the identity matrix
		numeric_matrix factor(arnoldiSize);
		exponential.copy_from(identityMatrix);

			factor = base;
			exponential.add(factor,1);

			for(int j=2;j<taylorSteps;j++) // Taylor series as matrix exponential
			{
				factor.mmult(factor,base);
				factor /=j;
				exponential.add (factor, 1);
			}
		return exponential;
	}


	numeric_matrix padeApproximation2 (numeric_matrix base)
	{
		int m = taylorSteps;

		double s = base.linfty_norm();

						if (s > 0.5)
						{

						s = trunc(log(s)/log(2))+2;
						if (s<0)
							{
							s=0;
							}
						 base *= (pow(2,(-s)));
						}



		numeric_vector c(m+2);
						c(1) = 1.0;
						for (int k = 1;k<=m;k++)
							{
							double a = k;
							double b = m;
							c(k+1) = c(k)*((b+1-a)/(a*(2*b+1-a)));
							}

						numeric_matrix I;
						I = makeIdentityMatrix();
						numeric_matrix aid(arnoldiSize);
						numeric_matrix E(arnoldiSize);
						numeric_matrix A2(arnoldiSize);
						numeric_matrix Q(arnoldiSize);
						numeric_matrix P(arnoldiSize);
						aid = base;

						base.mmult(A2,aid);
						Q = I;
						Q *= c(m+1);
						P = I;
						P *= c(m);
						int odd = 1;
						for (int k = m-1;k>=1;k--)
						{
						  if (odd == 1)
						    {
							  aid = Q;
							  aid.mmult(Q,A2);
							  aid = I;
							  aid *= c(k);
							  Q.add(aid,1);
						    }
						  else
						  {
						  aid = P;
						  aid.mmult(P,A2);
						  aid = I;
						  aid *= c(k);
						  P.add(aid,1);
						}
						  odd = 1-odd;
						}


						if (odd == 1)
						{
						  aid = Q;
						  aid.mmult(Q,base);
						  Q.add(P,-1);
						  E = I;
						  P.gauss_jordan();
						  aid = Q;
						  aid.mmult(Q,P);
						  Q *= 2;
						  E.add(Q,2);
						  E*=-1;
						}
						  else
						  {
							  aid = P;
							  aid.mmult(P,base);
							  Q.add(P,-1);
							  E = I;
							  P.gauss_jordan();
							  aid = Q;
							  aid.mmult(Q,P);
							  Q *= 2;
							  E.add(Q,2);

						  }




	for (int i=0;i<s;i++)
		{
		aid=E;
		aid.mmult(Q,E);
		E = Q;
		}
	return E;
	}


	numeric_matrix padeApproximation (numeric_matrix base)
	{
		numeric_matrix exponential(arnoldiSize);  // Builds the identity matrix
				numeric_matrix factor(arnoldiSize);
				numeric_matrix aid(arnoldiSize);
				numeric_matrix N(arnoldiSize);
				numeric_matrix D(arnoldiSize);
				int m=taylorSteps;

				double s = base.linfty_norm();

								if (s > 0.5)
								{

								s = trunc(log(s)/log(2))+2;
								if (s<0)
									{
									s=0;
									}
								 base *= (pow(2,(-s)));
								}

				double fact = 1.0;

				factor.copy_from(identityMatrix);
				N=factor;
				D=factor;





				for(int n=1;n<=m;n++) // Padé approximation as matrix exponential
							{
					double a = m;
					double b = n;

				fact = factorial(2*a-b)*factorial(a)/(factorial(2*a)*factorial(b)*factorial(a-b));

				factor.mmult(aid,base);
				factor = aid;

				aid *= fact;
				N.add(aid,1);

				if(n%2!=0)
				{
					aid *= (-1);
				}

				D.add(aid,1);

							}
				D.gauss_jordan();

				D.mmult(exponential,N);

				for (int i=1;i<=s;i++)
				{
					D=exponential;
					D.mmult(N,exponential);
					exponential = N;
				}

		return exponential;

	}

	void phiFunction(numeric_matrix H_m,numeric_matrix &phiOne, numeric_matrix &phiTwo) // Builds phiOne and phiTwo
	{
		numeric_matrix exponential(arnoldiSize);
		numeric_matrix inverse(arnoldiSize);
		phiOne = H_m;
		phiOne *= this->getTimeStepSize();

		exponential = padeApproximation(phiOne);
		exponential.add(makeIdentityMatrix(),(-1));

		phiOne.gauss_jordan();

		inverse = phiOne;

		exponential.mmult(phiOne,inverse);



		phiTwo = phiOne;
		phiTwo.add(makeIdentityMatrix(),(-1));

		exponential = phiTwo;
		exponential.mmult(phiTwo,inverse);

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
