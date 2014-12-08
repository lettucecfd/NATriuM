/**
 * @file ExponentialTimeIntegrator.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */
#include "ExponentialTimeIntegrator.h"
#include <iostream>
#include <cassert>
namespace natrium {

template<> ExponentialTimeIntegrator<distributed_sparse_matrix, distributed_vector>::ExponentialTimeIntegrator(
		double timeStepSize, size_t problemSize) :
		TimeIntegrator<distributed_sparse_matrix, distributed_vector>(
				timeStepSize), identityMatrix(makeIdentityMatrix()), H_m(makeMatrix(arnoldiSize,arnoldiSize)), H(makeMatrix(arnoldiSize+1,arnoldiSize))   {

}

template<> ExponentialTimeIntegrator<distributed_sparse_block_matrix,
		distributed_block_vector>::ExponentialTimeIntegrator(double timeStepSize,
		size_t problemSize, size_t numberOfBlocks) :
		TimeIntegrator<distributed_sparse_block_matrix, distributed_block_vector>(
				timeStepSize), identityMatrix(makeIdentityMatrix()), H_m(makeMatrix(arnoldiSize,arnoldiSize)), H(makeMatrix(arnoldiSize+1,arnoldiSize))  {
}



template<class MATRIX, class VECTOR> void ExponentialTimeIntegrator<MATRIX, VECTOR>::step(
		VECTOR& f, const MATRIX& systemMatrix, const VECTOR& systemVector) {

	// Test all dimensions and change, if necessary
	assert(systemMatrix.n() == systemMatrix.m());
	assert(f.size() == systemMatrix.n());

	if (w.size() != f.size()) {
			w.reinit(f.size(),true);
	}

	if (v_j.size() != f.size()) {
			v_j.reinit(f.size(),true);
	}

	if (v_i.size() != f.size()) {
			v_i.reinit(f.size(),true);
	}

	if (firstColumn.size() != arnoldiSize) {
			firstColumn.reinit(arnoldiSize,true);
		}

	if (m_f.size() != f.size()) {
				m_f.reinit(f.size(),true);
		}

	w = f;
	v_j = f;
	v_i = f;


numeric_matrix V(f.size(),arnoldiSize); // Orthonormal basis


double f_Norm = f.l2_norm();

for (int i=0;i<f.size();i++) // Arnoldi algorithm (first step)

{
	V.set(i,0,f(i) / f_Norm );

}

for (int j=0;j<arnoldiSize;j++)  // Arnoldi algorithm (second step)
	{


		for (int i=0;i<f.size();i++)
		{
			v_j(i) = V(i,j);

		}

		systemMatrix.vmult(w,v_j);
		w *= this->getTimeStepSize();

		for (int i=0;i<=j;i++)
		{
				for (int k=0;k<f.size();k++)
				{
					v_i(k) = V(k,i);
				}

			H(i,j) = w * v_i;
			v_i *= H(i,j);
			w -= v_i;

		}

		H(j+1,j) = w.l2_norm();
		if(H(j+1,j)!=0 && j<arnoldiSize-1)
		{
			for (int k=0;k<f.size();k++)
			{
				V(k,j+1) = w(k) / H(j+1,j);

			}
		}
	}

for (int i=0;i<arnoldiSize;i++) // Transform H to H_m (symmetric)
		for(int j=0;j<arnoldiSize;j++)
		{
			{
				H_m(i,j) = H(i,j);
			}
		}


H_m=taylorSeries(H_m); // Matrix exponential of H_m


for(int i=0;i<arnoldiSize;i++)
{
firstColumn(i)=H_m(i,0);
}


V.vmult(m_f,firstColumn);
m_f*=f_Norm;
f = m_f;



	}

template void ExponentialTimeIntegrator<distributed_sparse_matrix,
		distributed_vector>::step(distributed_vector& f,
		const distributed_sparse_matrix& systemMatrix, const distributed_vector& systemVector);
template void ExponentialTimeIntegrator<distributed_sparse_block_matrix,
		distributed_block_vector>::step(distributed_block_vector& f,
		const distributed_sparse_block_matrix& systemMatrix, const distributed_block_vector& systemVector);

} /* namespace natrium */
