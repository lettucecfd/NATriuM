/**
 * @file TimeIntegrator.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "TimeIntegrator.h"

namespace natrium {

template<class MATRIX, class VECTOR> TimeIntegrator<MATRIX, VECTOR>::TimeIntegrator(
		double timeStepSize) :
		m_timeStepSize(timeStepSize) {
}
template TimeIntegrator<distributed_sparse_matrix, distributed_vector>::TimeIntegrator(
		double timeStepSize);
template TimeIntegrator<distributed_sparse_block_matrix,
		distributed_block_vector>::TimeIntegrator(double timeStepSize);

/// destructor
//template <class MATRIX, class VECTOR> TimeIntegrator<MATRIX, VECTOR>::~TimeIntegrator(){
//}

} /* namespace natrium */
