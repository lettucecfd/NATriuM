/**
 * @file ExponentialTimeIntegrator.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */
#include "ExponentialTimeIntegrator.h"
#include <iostream>
#include <cassert>

#include "../utilities/Logging.h"

namespace natrium {

template<> ExponentialTimeIntegrator<distributed_sparse_matrix,
		distributed_vector>::ExponentialTimeIntegrator(double timeStepSize,
		size_t problemSize) :
		TimeIntegrator<distributed_sparse_matrix, distributed_vector>(
				timeStepSize), m_identityMatrix(makeIdentityMatrix()), m_Hm(
				makeMatrix(arnoldiSize, arnoldiSize)), m_H(
				makeMatrix(arnoldiSize + 2, arnoldiSize + 2)), m_phiOne(makeMatrix(arnoldiSize, arnoldiSize)),
				m_phiTwo(makeMatrix(arnoldiSize, arnoldiSize)), m_phiExtended(makeMatrix(arnoldiSize + 1,arnoldiSize + 1)) {

}

template<> ExponentialTimeIntegrator<distributed_sparse_block_matrix,
		distributed_block_vector>::ExponentialTimeIntegrator(
		double timeStepSize, size_t problemSize, size_t numberOfBlocks) :
		TimeIntegrator<distributed_sparse_block_matrix, distributed_block_vector>(
				timeStepSize), m_identityMatrix(makeIdentityMatrix()), m_Hm(
				makeMatrix(arnoldiSize, arnoldiSize)), m_H(
				makeMatrix(arnoldiSize + 2, arnoldiSize + 2)), m_phiOne(makeMatrix(arnoldiSize, arnoldiSize)),
				m_phiTwo(makeMatrix(arnoldiSize, arnoldiSize)), m_phiExtended(makeMatrix(arnoldiSize + 1,arnoldiSize + 1)) {
}


template<class MATRIX, class VECTOR> double ExponentialTimeIntegrator<MATRIX,
		VECTOR>::step(VECTOR& f, const MATRIX& systemMatrix,
		const VECTOR& systemVector, double t, double dt) {

	// Test all dimensions and change, if necessary
	assert(systemMatrix.n() == systemMatrix.m());
	assert(f.size() == systemMatrix.n());

	if ((0.0 != dt) and dt != this->getTimeStepSize()) {
		this->setTimeStepSize(dt);
		LOG(BASIC) << "Time step size set to " << dt << endl;
	} else if (0.0 == dt) {
		dt = this->getTimeStepSize();
	}


#ifdef WITH_TRILINOS
	if (m_w.size() != f.size()) {
		m_w.reinit(f);
	}

	if (m_vj.size() != f.size()) {
		m_vj.reinit(f);
	}

	if (m_vi.size() != f.size()) {
		m_vi.reinit(f);
	}
#else
	if (m_w.size() != f.size()) {
		m_w.reinit(f.size(), true);
	}

	if (m_vj.size() != f.size()) {
		m_vj.reinit(f.size(), true);
	}

	if (m_vi.size() != f.size()) {
		m_vi.reinit(f.size(), true);
	}
#endif

	if (m_firstColumn.size() != arnoldiSize + 1) {
		m_firstColumn.reinit(arnoldiSize + 1, true);
	}

	if (m_f.size() != f.size()) {
		m_f.reinit(f.size(), true);
	}



	m_w = f;
	m_vj = f;
	m_vi = f;

	for (int i = 0; i < arnoldiSize + 2; i++)
	{
		for (int j = 0; j < arnoldiSize + 2; j++)
		{
			m_H(i, j) = 0;
		}
	}

	numeric_matrix V(f.size(), arnoldiSize + 1); // Orthonormal basis V_(m+1)

	systemMatrix.vmult(m_w, f); 	// w = A*f + u
	m_w.add(systemVector);		// w = A*f + u

	double beta = m_w.l2_norm();

	for (int i = 0; i < f.size(); i++) 		// Arnoldi algorithm (first step)

	{
		V.set(i, 0, m_w(i) / beta);
	}

	for (int j = 0; j < arnoldiSize; j++)  	// Arnoldi algorithm (second step)
	{
		for (int i = 0; i < f.size(); i++)
		{
			m_vj(i) = V(i, j);
		}

		systemMatrix.vmult(m_w, m_vj);

		for (int i = 0; i <= j; i++)
		{
			for (int k = 0; k < f.size(); k++)
			{
				m_vi(k) = V(k, i);
			}

			m_H(i, j) = m_w * m_vi;
			m_vi *= m_H(i, j);
			m_w -= m_vi;

		}

		m_H(j + 1, j) = m_w.l2_norm();
		if (m_H(j + 1, j) != 0) // && j<arnoldiSize-1)
		{
			for (int k = 0; k < f.size(); k++)
				{
					V(k, j + 1) = m_w(k) / m_H(j + 1, j);
				}
		}
		m_H(arnoldiSize + 1, arnoldiSize) = 1;

	}

	for (int i = 0; i < arnoldiSize; i++)	// Transform H to Hm
	{
		for (int j = 0; j < arnoldiSize; j++)
		{
			m_Hm(i, j) = m_H(i, j);
		}
	}

/*	numeric_matrix phiOne(arnoldiSize);
	numeric_matrix phiTwo(arnoldiSize);
	numeric_matrix phiExtended;*/

	phiFunction(m_Hm, m_phiOne, m_phiTwo);

	for (int i = 0; i < arnoldiSize; i++)
	{
		for (int j = 0; j < arnoldiSize; j++)
		{
			m_phiExtended(i, j) = m_phiOne(i, j);
		}
		m_phiExtended(i, arnoldiSize) = 0;
	}

	for (int j = 0; j < arnoldiSize; j++)
	{
		m_phiExtended(arnoldiSize, j) = m_H(arnoldiSize, arnoldiSize - 1)
				* this->getTimeStepSize();
		m_phiExtended(arnoldiSize, j) *= m_phiTwo(arnoldiSize - 1, j);
	}

	m_phiExtended(arnoldiSize, arnoldiSize) = 1;

	for (int i = 0; i < arnoldiSize + 1; i++)
	{
		m_firstColumn(i) = m_phiExtended(i, 0);
	}
	V.vmult(m_f, m_firstColumn);

	m_f *= (dt * beta);

	for (int i = 0; i < f.size(); i++)

	{
		f(i) += m_f(i);
	}

	/*V.vmult(m_f,m_firstColumn);
	 m_f*=beta;
	 f = m_f;
	 */

	return t + dt;



}

template<class MATRIX, class VECTOR> double ExponentialTimeIntegrator<MATRIX,VECTOR>::factorial(int base)
{
	unsigned int factorial =1;
	for (int i=1;i<=base;i++)
	{
		factorial = factorial * i;
	}
	return factorial;
}

template<class MATRIX, class VECTOR> numeric_matrix ExponentialTimeIntegrator<MATRIX,VECTOR>:: makeIdentityMatrix()
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

template<class MATRIX, class VECTOR> numeric_matrix ExponentialTimeIntegrator<MATRIX,VECTOR>::makeMatrix(size_t m, size_t n)
	{
		numeric_matrix newMatrix(m,n);
		return newMatrix;
	}

// Calculates the matrix exponential of the given matrix;
template<class MATRIX, class VECTOR> numeric_matrix ExponentialTimeIntegrator<MATRIX,VECTOR>::taylorSeries(numeric_matrix base)
	{
		numeric_matrix exponential(arnoldiSize);  // Builds the identity matrix
		numeric_matrix factor(arnoldiSize);
		exponential.copy_from(m_identityMatrix);

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

template<class MATRIX, class VECTOR> numeric_matrix ExponentialTimeIntegrator<MATRIX,VECTOR>::padeApproximation (numeric_matrix base)
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

		factor.copy_from(m_identityMatrix);
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

template<class MATRIX, class VECTOR> void ExponentialTimeIntegrator<MATRIX,VECTOR>:: phiFunction(const numeric_matrix &Hm,numeric_matrix &phiOne, numeric_matrix &phiTwo) // Builds phiOne and phiTwo
	{
		numeric_matrix exponential(arnoldiSize);
		numeric_matrix inverse(arnoldiSize);
		phiOne = Hm;
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


template double ExponentialTimeIntegrator<distributed_sparse_matrix,
		distributed_vector>::step(distributed_vector& f,
		const distributed_sparse_matrix& systemMatrix,
		const distributed_vector& systemVector, double t, double dt);
template double ExponentialTimeIntegrator<distributed_sparse_block_matrix,
		distributed_block_vector>::step(distributed_block_vector& f,
		const distributed_sparse_block_matrix& systemMatrix,
		const distributed_block_vector& systemVector, double t, double dt);

template class ExponentialTimeIntegrator<distributed_sparse_block_matrix, distributed_block_vector>;
template class ExponentialTimeIntegrator<distributed_sparse_matrix, distributed_vector>;





} /* namespace natrium */
