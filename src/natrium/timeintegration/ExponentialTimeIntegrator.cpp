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
				timeStepSize) {

}

template<> ExponentialTimeIntegrator<distributed_sparse_block_matrix,
		distributed_block_vector>::ExponentialTimeIntegrator(double timeStepSize,
		size_t problemSize, size_t numberOfBlocks) :
		TimeIntegrator<distributed_sparse_block_matrix, distributed_block_vector>(
				timeStepSize) {
}



template<class MATRIX, class VECTOR> void ExponentialTimeIntegrator<MATRIX, VECTOR>::step(
		VECTOR& f, const MATRIX& systemMatrix, const VECTOR& systemVector) {

	// Test all dimensions and change, if necessary
	assert(systemMatrix.n() == systemMatrix.m());
	assert(f.size() == systemMatrix.n());


cout << "Dimension" << systemMatrix.n() << " f.size " << f.size()<< " " << f.l2_norm() <<  endl;

const int m = 5; // Factor of the arnoldi algorithm; sets the size of the orthonomal matrix V



numeric_matrix m_Adt(systemMatrix.m(),systemMatrix.n());
numeric_matrix V(f.size(),m); // Orthonormal basis
numeric_matrix H(m+1,m); // Hessenberg matrix
numeric_matrix H_m(m); // Hessenberg matrix (symmetric)
numeric_vector v_j(f.size());

m_Adt.copy_from(systemMatrix);
m_Adt *= this->getTimeStepSize();

double f_Norm = f.l2_norm();

for (int i=0;i<f.size();i++) // Arnoldi algorithm (first step)

{
	V(i,0) = f(i) / f_Norm ;

}

for (int j=0;j<m;j++)  // Arnoldi algorithm (second step)
	{


		for (int i=0;i<f.size();i++)
		{
			v_j(i) = V(i,j);

		}

		numeric_vector w(f.size());
		cout << "Test hier" << endl;
		m_Adt.vmult(w,v_j);
	//	w *=this->getTimeStepSize();

		for (int i=0;i<=j;i++)
		{
			numeric_vector v_i(f.size());
				for (int k=0;k<f.size();k++)
				{
					v_i(k) = V(k,i);
				}

			H(i,j) = w * v_i;
			v_i *= H(i,j);
			w -= v_i;

		}

		H(j+1,j) = w.l2_norm();
		if(H(j+1,j)!=0 && j<m-1)
		{
			for (int k=0;k<f.size();k++)
			{
				V(k,j+1) = w(k) / H(j+1,j);

			}
		}
	}

for (int i=0;i<m;i++)
		for(int j=0;j<m;j++)
		{
			{
				H_m(i,j) = H(i,j);
			}
		}


numeric_matrix exponential(m);  // Builds the identity matrix
	for (int i=0;i<m;i++)
		for(int j=0;j<m;j++)
		{
			{
				if(i==j)
					exponential(i,j)=1;
				else
					exponential(i,j)=0;
			}
		}



	numeric_matrix factor(H_m.m(),H_m.n());

	factor = H_m;
	exponential.add(factor,1);

	for(int j=2;j<29;j++) // Taylor series as matrix exponential
	{

		factor.mmult(factor,H_m);
		factor /=j;
		exponential.add (factor, 1);
	}


numeric_vector e(m);  // Vector e_1
for(int i=0;i<m;i++)
{
	if(i==0)
		e(i)=1;
	else
		e(i)=0;
}


numeric_vector f_f(f.size());  // auxiliary vector (to be removed)
numeric_vector aid(m);
exponential.vmult(aid,e);
//V.mmult(helpme,exponential);
V.vmult(f_f,aid);
f_f*=f_Norm;
f = f_f;

// = f_Norm * V * exponential*e;


	}

template void ExponentialTimeIntegrator<distributed_sparse_matrix,
		distributed_vector>::step(distributed_vector& f,
		const distributed_sparse_matrix& systemMatrix, const distributed_vector& systemVector);
template void ExponentialTimeIntegrator<distributed_sparse_block_matrix,
		distributed_block_vector>::step(distributed_block_vector& f,
		const distributed_sparse_block_matrix& systemMatrix, const distributed_block_vector& systemVector);

} /* namespace natrium */
