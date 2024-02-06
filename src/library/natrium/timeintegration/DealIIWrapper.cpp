/*
 * DealIIWrapper.cpp
 *
 *  Created on: Feb 5, 2015
 *      Author: kraemer
 */

#include "DealIIWrapper.h"

#include "deal.II/lac/solver_control.h"
#include "deal.II/lac/identity_matrix.h"

#include "deal.II/lac/solver_bicgstab.h"
#include "deal.II/lac/solver_cg.h"
#include "deal.II/lac/solver_gmres.h"
#include "deal.II/lac/solver_minres.h"
#include "deal.II/lac/solver_qmrs.h"
#include "deal.II/lac/solver_richardson.h"
#ifdef WITH_TRILINOS
#include "deal.II/lac/trilinos_precondition.h"
#include "../utilities/TrilinosBlockPreconditioner.h"
#endif

#include "deal.II/lac/precondition.h"

#include "../utilities/Logging.h"

namespace natrium {

/// TODO: refactor this class;  See class dealii::SolverSelector
/// Maybe with a new LinearSolver class
/// Unif

template<class MATRIX, class VECTOR>
VECTOR natrium::DealIIWrapper<MATRIX, VECTOR>::evaluateF(const double,
		const VECTOR& f) const {
	VECTOR result = *m_systemVector;

	// A*f(t) + b
	{
		TimerOutput::Scope timer_section(Timing::getTimer(), "vmult");
		m_systemMatrix->vmult_add(result, f);
	}

	return result;
}

template<>
distributed_vector natrium::DealIIWrapper<distributed_sparse_matrix,
		distributed_vector>::evaluateJInverse(const double, const double,
		const distributed_vector& f) const {

	distributed_vector result = f;

	dealii::SolverControl solver_control(m_iterations, m_tol * f.l2_norm(),
			false, false);	//* m_tmpSystemVector.l2_norm());

	switch (m_solver) {
	case 0: {

		dealii::SolverBicgstab<distributed_vector> bicgstab(solver_control);
		bicgstab.solve(*m_systemMatrix, result, f,
				dealii::TrilinosWrappers::PreconditionIdentity());
		break;
	}
	case 1: {
		dealii::SolverCG<distributed_vector> cg(solver_control);
		cg.solve(*m_systemMatrix, result, f,
				dealii::TrilinosWrappers::PreconditionIdentity());
		break;
	}
	case 2: {
		dealii::SolverFGMRES<distributed_vector> fgmres(solver_control);
		fgmres.solve(*m_systemMatrix, result, f,
				dealii::TrilinosWrappers::PreconditionIdentity());
		break;
	}
	case 3: {
		dealii::SolverGMRES<distributed_vector> gmres(solver_control);
		gmres.solve(*m_systemMatrix, result, f,
				dealii::TrilinosWrappers::PreconditionIdentity());
		break;
	}
	case 4: {
		dealii::SolverMinRes<distributed_vector> minres(solver_control);
		minres.solve(*m_systemMatrix, result, f,
				dealii::TrilinosWrappers::PreconditionIdentity());
		break;
	}
	case 5: {
		dealii::SolverQMRS<distributed_vector> qmrs(solver_control);
		qmrs.solve(*m_systemMatrix, result, f,
				dealii::TrilinosWrappers::PreconditionIdentity());
		break;
	}
	case 7: {
		dealii::SolverRichardson<distributed_vector> richard(solver_control);
		richard.solve(*m_systemMatrix, result, f,
				dealii::TrilinosWrappers::PreconditionIdentity());
		break;
	}

	default: {
		pout
				<< "Something went wrong with the selection of the Deal.II linear solver"
				<< endl;
	}
	}

	return result;
}

template<>
numeric_vector natrium::DealIIWrapper<sparse_matrix, numeric_vector>::evaluateJInverse(
		const double, const double, const numeric_vector& f) const {

	numeric_vector result = f;

	dealii::SolverControl solver_control(m_iterations, m_tol * f.l2_norm(),
			false, false);	//* m_tmpSystemVector.l2_norm());

	switch (m_solver) {
	case 0: {

		dealii::SolverBicgstab<numeric_vector> bicgstab(solver_control);
		bicgstab.solve(*m_systemMatrix, result, f,
				dealii::PreconditionIdentity());
		break;
	}
	case 1: {
		dealii::SolverCG<numeric_vector> cg(solver_control);
		cg.solve(*m_systemMatrix, result, f, dealii::PreconditionIdentity());
		break;
	}
	case 2: {
		dealii::SolverFGMRES<numeric_vector> fgmres(solver_control);
		fgmres.solve(*m_systemMatrix, result, f,
				dealii::PreconditionIdentity());
		break;
	}
	case 3: {
		dealii::SolverGMRES<numeric_vector> gmres(solver_control);
		gmres.solve(*m_systemMatrix, result, f, dealii::PreconditionIdentity());
		break;
	}
	case 4: {
		dealii::SolverMinRes<numeric_vector> minres(solver_control);
		minres.solve(*m_systemMatrix, result, f,
				dealii::PreconditionIdentity());
		break;
	}
	case 5: {
		dealii::SolverQMRS<numeric_vector> qmrs(solver_control);
		qmrs.solve(*m_systemMatrix, result, f, dealii::PreconditionIdentity());
		break;
	}
	case 7: {
		dealii::SolverRichardson<numeric_vector> richard(solver_control);
		richard.solve(*m_systemMatrix, result, f,
				dealii::PreconditionIdentity());
		break;
	}
	}	//,	           preconditioner);
	return result;
}

template<>
distributed_block_vector natrium::DealIIWrapper<distributed_sparse_block_matrix,
		distributed_block_vector>::evaluateJInverse(const double, const double,
		const distributed_block_vector& f) const {

	distributed_block_vector result(f);

	dealii::SolverControl solver_control(m_iterations, m_tol * f.l2_norm(),
			false, false);	//* m_tmpSystemVector.l2_norm());

	switch (m_solver) {
	case 0: {

		dealii::SolverBicgstab<distributed_block_vector> bicgstab(
				solver_control);
		bicgstab.solve(*m_systemMatrix, result, f,
				TrilinosBlockPreconditioner());
		break;
	}
	case 1: {
		dealii::SolverCG<distributed_block_vector> cg(solver_control);
		cg.solve(*m_systemMatrix, result, f, TrilinosBlockPreconditioner());
		break;
	}
	case 2: {
		dealii::SolverFGMRES<distributed_block_vector> fgmres(solver_control);
		fgmres.solve(*m_systemMatrix, result, f, TrilinosBlockPreconditioner());
		break;
	}
	case 3: {
		dealii::SolverGMRES<distributed_block_vector> gmres(solver_control);
		gmres.solve(*m_systemMatrix, result, f, TrilinosBlockPreconditioner());
		break;
	}
	case 4: {
		dealii::SolverMinRes<distributed_block_vector> minres(solver_control);
		minres.solve(*m_systemMatrix, result, f, TrilinosBlockPreconditioner());
		break;
	}
	case 5: {
		dealii::SolverQMRS<distributed_block_vector> qmrs(solver_control);
		qmrs.solve(*m_systemMatrix, result, f, TrilinosBlockPreconditioner());
		break;
	}
	case 7: {
		dealii::SolverRichardson<distributed_block_vector> richard(
				solver_control);
		richard.solve(*m_systemMatrix, result, f,
				TrilinosBlockPreconditioner());
		break;
	}

	default: {
		pout
				<< "Something went wrong with the selection of the Deal.II linear solver"
				<< endl;
	}
	}

	return result;
}

template<>
block_vector natrium::DealIIWrapper<sparse_block_matrix, block_vector>::evaluateJInverse(
		const double, const double, const block_vector& f) const {

	block_vector result(f);

	dealii::SolverControl solver_control(m_iterations, m_tol * f.l2_norm(),
			false, false);	//* m_tmpSystemVector.l2_norm());

	switch (m_solver) {
	case 0: {

		dealii::SolverBicgstab<block_vector> bicgstab(solver_control);
		bicgstab.solve(*m_systemMatrix, result, f,
				dealii::PreconditionIdentity());
		break;
	}
	case 1: {
		dealii::SolverCG<block_vector> cg(solver_control);
		cg.solve(*m_systemMatrix, result, f, dealii::PreconditionIdentity());
		break;
	}
	case 2: {
		dealii::SolverFGMRES<block_vector> fgmres(solver_control);
		fgmres.solve(*m_systemMatrix, result, f,
				dealii::PreconditionIdentity());
		break;
	}
	case 3: {
		dealii::SolverGMRES<block_vector> gmres(solver_control);
		gmres.solve(*m_systemMatrix, result, f, dealii::PreconditionIdentity());
		break;
	}
	case 4: {
		dealii::SolverMinRes<block_vector> minres(solver_control);
		minres.solve(*m_systemMatrix, result, f,
				dealii::PreconditionIdentity());
		break;
	}
	case 5: {
		dealii::SolverQMRS<block_vector> qmrs(solver_control);
		qmrs.solve(*m_systemMatrix, result, f, dealii::PreconditionIdentity());
		break;
	}
	case 7: {
		dealii::SolverRichardson<block_vector> richard(solver_control);
		richard.solve(*m_systemMatrix, result, f,
				dealii::PreconditionIdentity());
		break;
	}
	}	//,	           preconditioner);
	return result;
}

template<>
distributed_vector natrium::DealIIWrapper<distributed_sparse_matrix,
		distributed_vector>::evaluateIdMinusTauJInverse(const double,
		const double tau, const distributed_vector& f) {

	distributed_vector result = f;

	// A
	m_tmpMatrix.copy_from(*m_systemMatrix);
	// I-theta*A
	m_tmpMatrix *= (-tau);
	for (size_t i = 0; i < m_tmpMatrix.n(); i++) {
		m_tmpMatrix.add(i, i, 1.0);
	}

	dealii::SolverControl solver_control(m_iterations, m_tol * f.l2_norm(),
			false, false);	//* m_tmpSystemVector.l2_norm());

	switch (m_solver) {
	case 0: {

		dealii::SolverBicgstab<distributed_vector> bicgstab(solver_control);
		bicgstab.solve(m_tmpMatrix, result, f,
				dealii::TrilinosWrappers::PreconditionIdentity());
		break;
	}
	case 1: {
		dealii::SolverCG<distributed_vector> cg(solver_control);
		cg.solve(m_tmpMatrix, result, f,
				dealii::TrilinosWrappers::PreconditionIdentity());
		break;
	}
	case 2: {
		dealii::SolverFGMRES<distributed_vector> fgmres(solver_control);
		fgmres.solve(m_tmpMatrix, result, f,
				dealii::TrilinosWrappers::PreconditionIdentity());
		break;
	}
	case 3: {
		dealii::SolverGMRES<distributed_vector> gmres(solver_control);
		gmres.solve(m_tmpMatrix, result, f,
				dealii::TrilinosWrappers::PreconditionIdentity());
		break;
	}
	case 4: {
		dealii::SolverMinRes<distributed_vector> minres(solver_control);
		minres.solve(m_tmpMatrix, result, f,
				dealii::TrilinosWrappers::PreconditionIdentity());
		break;
	}
	case 5: {
		dealii::SolverQMRS<distributed_vector> qmrs(solver_control);
		qmrs.solve(m_tmpMatrix, result, f,
				dealii::TrilinosWrappers::PreconditionIdentity());
		break;
	}
	case 7: {
		dealii::SolverRichardson<distributed_vector> richard(solver_control);
		richard.solve(m_tmpMatrix, result, f,
				dealii::TrilinosWrappers::PreconditionIdentity());
		break;
	}

	default: {
		pout
				<< "Something went wrong with the selection of the Deal.II linear solver"
				<< endl;
	}
	}

	return result;
}

template<>
numeric_vector natrium::DealIIWrapper<sparse_matrix, numeric_vector>::evaluateIdMinusTauJInverse(
		const double, const double tau, const numeric_vector& f) {

	numeric_vector result = f;

	if ((m_tmpMatrix.empty()) or
	// the next check should give true if the sparsity patterns are equal
	// and false, else. n_nonzero_elements returns the number of entries
	// in the sparsity pattern, not the actual number of nonzero entries
			(m_tmpMatrix.n_nonzero_elements()
					!= m_systemMatrix->n_nonzero_elements())) {
		m_tmpMatrix.reinit(m_systemMatrix->get_sparsity_pattern());
	}

	// A
	m_tmpMatrix.copy_from(*m_systemMatrix);
	// I-theta*A
	m_tmpMatrix *= (-tau);
	for (size_t i = 0; i < m_tmpMatrix.n(); i++) {
		m_tmpMatrix.add(i, i, 1.0);
	}

	dealii::SolverControl solver_control(m_iterations, m_tol * f.l2_norm(),
			false, false);	//* m_tmpSystemVector.l2_norm());

	switch (m_solver) {
	case 0: {

		dealii::SolverBicgstab<numeric_vector> bicgstab(solver_control);
		bicgstab.solve(m_tmpMatrix, result, f, dealii::PreconditionIdentity());
		break;
	}
	case 1: {
		dealii::SolverCG<numeric_vector> cg(solver_control);
		cg.solve(m_tmpMatrix, result, f, dealii::PreconditionIdentity());
		break;
	}
	case 2: {
		dealii::SolverFGMRES<numeric_vector> fgmres(solver_control);
		fgmres.solve(m_tmpMatrix, result, f, dealii::PreconditionIdentity());
		break;
	}
	case 3: {
		dealii::SolverGMRES<numeric_vector> gmres(solver_control);
		gmres.solve(m_tmpMatrix, result, f, dealii::PreconditionIdentity());
		break;
	}
	case 4: {
		dealii::SolverMinRes<numeric_vector> minres(solver_control);
		minres.solve(m_tmpMatrix, result, f, dealii::PreconditionIdentity());
		break;
	}
	case 5: {
		dealii::SolverQMRS<numeric_vector> qmrs(solver_control);
		qmrs.solve(m_tmpMatrix, result, f, dealii::PreconditionIdentity());
		break;
	}
	case 7: {
		dealii::SolverRichardson<numeric_vector> richard(solver_control);
		richard.solve(m_tmpMatrix, result, f, dealii::PreconditionIdentity());
		break;
	}
	}	//,	           preconditioner);

	return result;
}

template<>
distributed_block_vector natrium::DealIIWrapper<distributed_sparse_block_matrix,
		distributed_block_vector>::evaluateIdMinusTauJInverse(const double,
		const double tau, const distributed_block_vector& f) {
	// Test all dimensions and change, if necessary
	assert(m_systemMatrix->n() == m_systemMatrix->m());
	assert(f.size() == m_systemMatrix->n());

	distributed_block_vector result = f;

	// check equality of sparsity patterns
	if (m_tmpMatrix.memory_consumption()
			!= m_systemMatrix->memory_consumption()) {
		size_t n_blocks = m_systemMatrix->n_block_rows();
		assert(
				m_systemMatrix->n_block_cols()
						== m_systemMatrix->n_block_rows());
		m_tmpMatrix.reinit(n_blocks, n_blocks);
		for (size_t I = 0; I < n_blocks; I++) {
			for (size_t J = 0; J < n_blocks; J++) {
				m_tmpMatrix.block(I, J).reinit(m_systemMatrix->block(I, J));
			}
		}
	}

	// A*f(t)
	m_tmpMatrix.copy_from(*m_systemMatrix);
	size_t n_blocks = m_systemMatrix->n_block_rows();
	for (size_t I = 0; I < n_blocks; I++) {
		for (size_t J = 0; J < n_blocks; J++) {
			assert(
					m_tmpMatrix.block(I, J).l1_norm()
							== m_systemMatrix->block(I, J).l1_norm());
		}
	}
	// I-tau*A
	m_tmpMatrix *= (-tau);
	const dealii::IndexSet& locally_owned_domain_indices = m_tmpMatrix.block(0,
			0).locally_owned_domain_indices();
	//for all degrees of freedom on current processor
	dealii::IndexSet::ElementIterator it(locally_owned_domain_indices.begin());
	dealii::IndexSet::ElementIterator end(locally_owned_domain_indices.end());
	for (size_t I = 0; I < m_tmpMatrix.n_block_cols(); I++) {
		for (it = locally_owned_domain_indices.begin(); it != end; it++) {
			size_t i = *it;
			m_tmpMatrix.block(I, I).add(i, i, 1.0);
		}
	}

	dealii::SolverControl solver_control(m_iterations, (m_tol * f.l2_norm()),
			false, false);	//* m_tmpSystemVector.l2_norm());

	switch (m_solver) {
	case 0: {
		dealii::SolverBicgstab<distributed_block_vector> bicgstab(
				solver_control);
		bicgstab.solve(m_tmpMatrix, result, f, dealii::PreconditionIdentity());
		break;
	}
	case 1: {
		dealii::SolverCG<distributed_block_vector> cg(solver_control);
		cg.solve(m_tmpMatrix, result, f, dealii::PreconditionIdentity());
		break;
	}
	case 2: {
		dealii::SolverFGMRES<distributed_block_vector> fgmres(solver_control);
		fgmres.solve(m_tmpMatrix, result, f, dealii::PreconditionIdentity());
		break;
	}
	case 3: {
		dealii::SolverGMRES<distributed_block_vector> gmres(solver_control);
		gmres.solve(m_tmpMatrix, result, f, dealii::PreconditionIdentity());
		break;
	}
	case 4: {
		dealii::SolverMinRes<distributed_block_vector> minres(solver_control);
		minres.solve(m_tmpMatrix, result, f, dealii::PreconditionIdentity());
		break;
	}
	case 5: {
		dealii::SolverQMRS<distributed_block_vector> qmrs(solver_control);
		qmrs.solve(m_tmpMatrix, result, f, dealii::PreconditionIdentity());
		break;
	}
	case 7: {
		dealii::SolverRichardson<distributed_block_vector> richard(
				solver_control);
		richard.solve(m_tmpMatrix, result, f, dealii::PreconditionIdentity());
		break;
	}

	default: {
		pout
				<< "Something went wrong with the selection of the Deal.II linear solver"
				<< endl;
	}
	}

	return result;

}

template<>
block_vector natrium::DealIIWrapper<sparse_block_matrix, block_vector>::evaluateIdMinusTauJInverse(
		const double, const double tau, const block_vector& f) {
	// Test all dimensions and change, if necessary
	assert(m_systemMatrix->n() == m_systemMatrix->m());
	assert(f.size() == m_systemMatrix->n());

	block_vector result = f;

	if ((m_tmpMatrix.empty()) or
	// the next check should give true if the sparsity patterns are equal
	// and false, else. n_nonzero_elements returns the number of entries
	// in the sparsity pattern, not the actual number of nonzero entries
			(m_tmpMatrix.n_nonzero_elements()
					!= m_systemMatrix->n_nonzero_elements())) {
		m_tmpMatrix.reinit(m_systemMatrix->get_sparsity_pattern());
	}

	// A*f(t)
	m_tmpMatrix.copy_from(*m_systemMatrix);
	size_t n_blocks = m_systemMatrix->n_block_rows();
	for (size_t I = 0; I < n_blocks; I++) {
		for (size_t J = 0; J < n_blocks; J++) {
			assert(
					m_tmpMatrix.block(I, J).l1_norm()
							== m_systemMatrix->block(I, J).l1_norm());
		}
	}
	// I-tau*A
	m_tmpMatrix *= (-tau);
	for (size_t I = 0; I < m_tmpMatrix.n_block_cols(); I++) {
		for (size_t i = 0; i < m_tmpMatrix.block(I, I).n(); i++) {
			m_tmpMatrix.block(I, I).add(i, i, 1.0);
		}
	}

	dealii::SolverControl solver_control(m_iterations, (m_tol * f.l2_norm()),
			false, false);	//* m_tmpSystemVector.l2_norm());

	switch (m_solver) {
	case 0: {

		dealii::SolverBicgstab<block_vector> bicgstab(solver_control);
		bicgstab.solve(m_tmpMatrix, result, f, dealii::PreconditionIdentity());

		break;
	}
	case 1: {
		dealii::SolverCG<block_vector> cg(solver_control);
		cg.solve(m_tmpMatrix, result, f, dealii::PreconditionIdentity());

		break;
	}
	case 2: {
		dealii::SolverFGMRES<block_vector> fgmres(solver_control);
		fgmres.solve(m_tmpMatrix, result, f, dealii::PreconditionIdentity());
		break;
	}
	case 3: {
		dealii::SolverGMRES<block_vector> gmres(solver_control);
		gmres.solve(m_tmpMatrix, result, f, dealii::PreconditionIdentity());
		break;
	}
	case 4: {
		dealii::SolverMinRes<block_vector> minres(solver_control);
		minres.solve(m_tmpMatrix, result, f, dealii::PreconditionIdentity());
		break;
	}
	case 5: {
		dealii::SolverQMRS<block_vector> qmrs(solver_control);
		qmrs.solve(m_tmpMatrix, result, f, dealii::PreconditionIdentity());
		break;
	}
	case 7: {
		dealii::SolverRichardson<block_vector> richard(solver_control);
		richard.solve(m_tmpMatrix, result, f, dealii::PreconditionIdentity());
		break;
	}
	}	//,	           preconditioner);

	return result;

}

template<class MATRIX, class VECTOR>
double natrium::DealIIWrapper<MATRIX, VECTOR>::step(VECTOR& vector,
		const MATRIX& systemMatrix, VECTOR& systemVector, double t,
		double dt) {
    (void) vector;
    (void) systemMatrix;
    (void) systemVector;
    (void) t;
    (void) dt;
    /*
TODO BRUTE FORCE KILLED SEDG HERE DOMINIK WILDE 2020

	if ((0.0 != dt) and dt != this->getTimeStepSize()) {
		this->setTimeStepSize(dt);
		if (m_dealIIRKEmbedded) {
			LOG(BASIC) << "Time step size set to " << dt << endl;
		} else {
			LOG(DETAILED) << "Time step size set to " << dt << endl;
		}
	} else if (0.0 == dt) {
		dt = this->getTimeStepSize();
	}

	m_systemMatrix = &systemMatrix;
	m_systemVector = &systemVector;
	return m_dealIIRKStepper->evolve_one_time_step(
			dealii::std_cxx11::bind(&DealIIWrapper<MATRIX, VECTOR>::evaluateF,
					this, dealii::std_cxx11::_1, dealii::std_cxx11::_2),
			dealii::std_cxx11::bind(
					&DealIIWrapper<MATRIX, VECTOR>::evaluateIdMinusTauJInverse,
					this, dealii::std_cxx11::_1, dealii::std_cxx11::_2,
					dealii::std_cxx11::_3), t, dt, vector);

					*/
    return 0;
}

template<class MATRIX, class VECTOR>
natrium::DealIIWrapper<MATRIX, VECTOR>::DealIIWrapper(const double timeStepSize,
		const DealIntegratorName rkScheme, const DealSolverName linearSolver) :
		TimeIntegrator<MATRIX, VECTOR>(timeStepSize), m_systemMatrix(NULL), m_systemVector(
		NULL), m_solver() {
	m_solver = linearSolver;
	assert(rkScheme < 12);
	assert(NONE != rkScheme);
	dealii::TimeStepping::runge_kutta_method rk =
			static_cast<dealii::TimeStepping::runge_kutta_method>((int) rkScheme);
	if (rkScheme < 3) {
		m_dealIIRKStepper = boost::make_shared<
				dealii::TimeStepping::ExplicitRungeKutta<VECTOR> >(rk);
	} else if (rkScheme < 7) {
		m_dealIIRKStepper = boost::make_shared<
				dealii::TimeStepping::ImplicitRungeKutta<VECTOR> >(rk);
	} else if (rkScheme < 12) {
		m_dealIIRKEmbedded = boost::make_shared<
				dealii::TimeStepping::EmbeddedExplicitRungeKutta<VECTOR> >(rk);
		m_dealIIRKStepper = m_dealIIRKEmbedded;
	};
}

template<class MATRIX, class VECTOR>
natrium::DealIIWrapper<MATRIX, VECTOR>::DealIIWrapper(const double timeStepSize,
		const DealIntegratorName rkScheme, const DealSolverName linearSolver,
		double coarsen_param, double refine_param, double min_delta,
		double max_delta, double refine_tol, double coarsen_tol) :
		TimeIntegrator<MATRIX, VECTOR>(timeStepSize), m_systemMatrix(NULL), m_systemVector(
		NULL), m_solver() {
	m_solver = linearSolver;
	assert(rkScheme < 12);
	assert(NONE != rkScheme);
	dealii::TimeStepping::runge_kutta_method rk =
			static_cast<dealii::TimeStepping::runge_kutta_method>((int) rkScheme);
	m_dealIIRKEmbedded = boost::make_shared<
			dealii::TimeStepping::EmbeddedExplicitRungeKutta<VECTOR> >(rk,
			coarsen_param, refine_param, min_delta, max_delta, refine_tol,
			coarsen_tol);
	m_dealIIRKStepper = m_dealIIRKEmbedded;

}

/// explicit instantiation
template class DealIIWrapper<distributed_sparse_matrix, distributed_vector> ;
template class DealIIWrapper<distributed_sparse_block_matrix,
		distributed_block_vector> ;
template class DealIIWrapper<sparse_matrix, numeric_vector> ;
template class DealIIWrapper<sparse_block_matrix, block_vector> ;
} /* namespace natrium */
