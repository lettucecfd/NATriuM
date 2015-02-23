/*
 * DealIIWrapper.cpp
 *
 *  Created on: Feb 5, 2015
 *      Author: kraemer
 */

#include <timeintegration/DealIIWrapper.h>

#include "deal.II/lac/solver_control.h"
#include "deal.II/lac/identity_matrix.h"
#include "deal.II/lac/solver_bicgstab.h"

#ifdef WITH_TRILINOS
#include "deal.II/lac/trilinos_precondition.h"
#include "../utilities/TrilinosBlockPreconditioner.h"
#else
#include "deal.II/lac/precondition.h"
#endif

namespace natrium {

template<class MATRIX, class VECTOR>
VECTOR natrium::DealIIWrapper<MATRIX, VECTOR>::evaluateF(const double t,
		const VECTOR& f) const {
	VECTOR result = *m_systemVector;

	// A*f(t) + b
	m_systemMatrix->vmult_add(result, f);

	return result;
}

template<>
distributed_vector natrium::DealIIWrapper<distributed_sparse_matrix,
		distributed_vector>::evaluateJInverse(const double t, const double tau,
		const distributed_vector& f) const {

	distributed_vector result = f;

	dealii::SolverControl solver_control(1000, 1e-8, false, false);	//* m_tmpSystemVector.l2_norm());
	dealii::SolverBicgstab<distributed_vector> bicgstab(solver_control);
#ifdef WITH_TRILINOS
	bicgstab.solve(*m_systemMatrix, result, f,
			dealii::TrilinosWrappers::PreconditionIdentity());
#else
	bicgstab.solve(*m_systemMatrix, result, f,
			dealii::PreconditionIdentity());	//,	           preconditioner);
#endif
	return result;
}

template<>
distributed_block_vector natrium::DealIIWrapper<distributed_sparse_block_matrix,
		distributed_block_vector>::evaluateJInverse(const double t,
		const double tau, const distributed_block_vector& f) const {

	distributed_block_vector result = f;

	dealii::SolverControl solver_control(1000, 1e-8, false, false);	//* m_tmpSystemVector.l2_norm());
	dealii::SolverBicgstab<distributed_block_vector> bicgstab(solver_control);
#ifdef WITH_TRILINOS
	bicgstab.solve(*m_systemMatrix, result, f, TrilinosBlockPreconditioner());
#else
	bicgstab.solve(*m_systemMatrix, result, f,
			dealii::PreconditionIdentity());	//,	           preconditioner);
#endif
	return result;
}

template<>
distributed_vector natrium::DealIIWrapper<distributed_sparse_matrix,
		distributed_vector>::evaluateIdMinusTauJInverse(const double t,
		const double tau, const distributed_vector& f) const {

	distributed_vector result = f;
	distributed_sparse_matrix tmpMatrix;
	// A
	tmpMatrix.copy_from(*m_systemMatrix);
	// I-theta*A
	tmpMatrix *= (-tau);
	for (size_t i = 0; i < tmpMatrix.n(); i++) {
		tmpMatrix.add(i, i, 1.0);
	}

	dealii::SolverControl solver_control(1000, 1e-8, false, false);	//* m_tmpSystemVector.l2_norm());
	dealii::SolverBicgstab<distributed_vector> bicgstab(solver_control);
#ifdef WITH_TRILINOS
	bicgstab.solve(tmpMatrix, result, f,
			dealii::TrilinosWrappers::PreconditionIdentity());
#else
	bicgstab.solve(tmpMatrix, result, f,
			dealii::PreconditionIdentity());	//,	           preconditioner);
#endif

	return result;
}

template<>
distributed_block_vector natrium::DealIIWrapper<distributed_sparse_block_matrix,
		distributed_block_vector>::evaluateIdMinusTauJInverse(const double t,
		const double tau, const distributed_block_vector& f) const {

	distributed_block_vector result = f;
	distributed_sparse_block_matrix tmpMatrix;
	// A
	tmpMatrix.copy_from(*m_systemMatrix);
	// I-theta*A
	tmpMatrix *= (-tau);
	for (size_t I = 0; I < tmpMatrix.n_block_cols(); I++) {
		for (size_t i = 0; i < tmpMatrix.block(I, I).n(); i++) {
			tmpMatrix.block(I, I).add(i, i, 1.0);
		}
	}

	dealii::SolverControl solver_control(1000, 1e-8, false, false);	//* m_tmpSystemVector.l2_norm());
	dealii::SolverBicgstab<distributed_block_vector> bicgstab(solver_control);
#ifdef WITH_TRILINOS
	bicgstab.solve(tmpMatrix, result, f, TrilinosBlockPreconditioner());
#else
	bicgstab.solve(tmpMatrix, result, f,
			dealii::PreconditionIdentity());	//,	           preconditioner);
#endif

	return result;
}

template<class MATRIX, class VECTOR>
void natrium::DealIIWrapper<MATRIX, VECTOR>::step(VECTOR& vector,
		const MATRIX& systemMatrix, const VECTOR& systemVector) {
	m_systemMatrix = &systemMatrix;
	m_systemVector = &systemVector;
	double t = 0;
	double dt = TimeIntegrator<MATRIX, VECTOR>::getTimeStepSize();
	m_dealIIRKStepper->evolve_one_time_step(
			dealii::std_cxx11::bind(&DealIIWrapper<MATRIX, VECTOR>::evaluateF,
					this, dealii::std_cxx11::_1, dealii::std_cxx11::_2),
			dealii::std_cxx11::bind(
					&DealIIWrapper<MATRIX, VECTOR>::evaluateIdMinusTauJInverse,
					this, dealii::std_cxx11::_1, dealii::std_cxx11::_2,
					dealii::std_cxx11::_3), t, dt, vector);
}

template<class MATRIX, class VECTOR>
natrium::DealIIWrapper<MATRIX, VECTOR>::DealIIWrapper(const double timeStepSize) :
		TimeIntegrator<MATRIX, VECTOR>(timeStepSize) {
	m_dealIIRKStepper = make_shared<
			dealii::TimeStepping::ExplicitRungeKutta<VECTOR> >(
			dealii::TimeStepping::RK_CLASSIC_FOURTH_ORDER);
}

/// explicit instantiation
template class DealIIWrapper<distributed_sparse_matrix, distributed_vector> ;
template class DealIIWrapper<distributed_sparse_block_matrix,
		distributed_block_vector> ;

} /* namespace natrium */
