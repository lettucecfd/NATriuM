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
#include "deal.II/lac/solver_relaxation.h"
#include "deal.II/lac/solver_richardson.h"


#ifdef WITH_TRILINOS
#include "deal.II/lac/trilinos_precondition.h"
#include "../utilities/TrilinosBlockPreconditioner.h"
#endif

#include "deal.II/lac/precondition.h"

#include "../utilities/Logging.h"

namespace natrium {


template<class MATRIX, class VECTOR>
VECTOR natrium::DealIIWrapper<MATRIX, VECTOR>::evaluateF(const double t,
        const VECTOR& f) const {
    VECTOR result = *m_systemVector;

    // A*f(t) + b
    m_systemMatrix->vmult_add(result, f);

    return result;
}

//template<class MATRIX, class VECTOR>
template<class MATRIX, class VECTOR, class PRECONDITION()>
VECTOR natrium::DealIIWrapper<MATRIX, VECTOR>::solvers(const VECTOR& f, bool isTmpMatrix)
        {

	cout << m_tmpMatrix.m() << " " << m_tmpMatrix.n() << " " << endl;
	cout << systemMatrix.memory_consumption() << endl;
	MATRIX testMatrix;
	if (isTmpMatrix)
	{
		//m_tmpMatrix.reinit();
		testMatrix.copy_from(m_tmpMatrix);
	}
	else
	{
		testMatrix.copy_from(*m_systemMatrix);
	}

	VECTOR result;
        dealii::SolverControl solver_control(1000, 1e-6*f.l2_norm(), false, false);	//* m_tmpSystemVector.l2_norm());

        switch(m_solver){

                case 0:
                {
                	cout << "point 10" << endl;
                    dealii::SolverBicgstab<VECTOR> bicgstab(solver_control);
                    bicgstab.solve(testMatrix, result, f, PRECONDITION());
                    cout << "BIC works" ;
                    break;
                }
                case 1:
                {
                    dealii::SolverCG<VECTOR> cg(solver_control);
                    cg.solve(testMatrix, result, f, PRECONDITION());
                    cout << "CG works" ;
                    break;
                }
                case 2:
                {
                    dealii::SolverFGMRES<VECTOR> fgmres(solver_control);
                    fgmres.solve(testMatrix, result, f,PRECONDITION());
                    break;
                }
                case 3:
                {
                    dealii::SolverGMRES<VECTOR> gmres(solver_control);
                    gmres.solve(testMatrix, result, f, PRECONDITION());
                    break;
                }
                case 4:
                {
                    dealii::SolverMinRes<VECTOR> minres(solver_control);
                    minres.solve(testMatrix, result, f, PRECONDITION());
                    break;
                }
                case 5:
                {
                    dealii::SolverQMRS<VECTOR> qmrs(solver_control);
                    qmrs.solve(testMatrix, result, f, PRECONDITION());
                    break;
                }

                case 7:
                {
                    dealii::SolverRichardson<VECTOR> richard(solver_control);
                    richard.solve(testMatrix, result, f, PRECONDITION());
                    break;
                }

                default:
                {
                    cout << "Something went wrong with the selection of the Deal.II linear solver" << endl;
                }
                }
        return result;



        }




template<>
distributed_vector natrium::DealIIWrapper<distributed_sparse_matrix,
        distributed_vector>::evaluateJInverse(const double t, const double tau,
        const distributed_vector& f) const {

    distributed_vector result = f;

#ifdef WITH_TRILINOS
    cout << "Point 1" << endl;
    result = solvers<dealii::TrilinosWrappers::PreconditionIdentity>(f, 0);

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


#ifdef WITH_TRILINOS
    cout << "Point 2" << endl;

    result = natrium::DealIIWrapper<distributed_sparse_block_matrix,
            distributed_block_vector>::solvers<distributed_sparse_block_matrix,distributed_block_vector,TrilinosBlockPreconditioner>(f, 0);

#else
    bicgstab.solve(*m_systemMatrix, result, f,
            dealii::PreconditionIdentity());	//,	           preconditioner);
#endif
    return result;
}

template<>
distributed_vector natrium::DealIIWrapper<distributed_sparse_matrix,
        distributed_vector>::evaluateIdMinusTauJInverse(const double t,
        const double tau, const distributed_vector& f) {

    distributed_vector result = f;

#ifndef WITH_TRILINOS
    if ((m_tmpMatrix.empty()) or
            // the next check should give true if the sparsity patterns are equal
            // and false, else. n_nonzero_elements returns the number of entries
            // in the sparsity pattern, not the actual number of nonzero entries
            (m_tmpMatrix.n_nonzero_elements()
                    != m_systemMatrix->n_nonzero_elements())) {
        m_tmpMatrix.reinit(m_systemMatrix->get_sparsity_pattern());
    }
#endif

    // A
    m_tmpMatrix.copy_from(*m_systemMatrix);
    // I-theta*A
    m_tmpMatrix *= (-tau);
    for (size_t i = 0; i < m_tmpMatrix.n(); i++) {
        m_tmpMatrix.add(i, i, 1.0);
    }


#ifdef WITH_TRILINOS
    result = solvers<dealii::TrilinosWrappers::PreconditionIdentity>(f, 1);

    #else
    bicgstab.solve(m_tmpMatrix, result, f,
            dealii::PreconditionIdentity());	//,	           preconditioner);
    #endif

    return result;
}

template<>
distributed_block_vector natrium::DealIIWrapper<distributed_sparse_block_matrix,
        distributed_block_vector>::evaluateIdMinusTauJInverse(const double t,
        const double tau, const distributed_block_vector& f) {
    // Test all dimensions and change, if necessary

	cout << m_tmpMatrix.memory_consumption() << endl;
	assert(m_systemMatrix->n() == m_systemMatrix->m());
    assert(f.size() == m_systemMatrix->n());

    distributed_block_vector result = f;
    cout << "Point 4" << endl;
#ifdef WITH_TRILINOS
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
#else
    if ((m_tmpMatrix.empty()) or
            // the next check should give true if the sparsity patterns are equal
            // and false, else. n_nonzero_elements returns the number of entries
            // in the sparsity pattern, not the actual number of nonzero entries
            (m_tmpMatrix.n_nonzero_elements()
                    != m_systemMatrix->n_nonzero_elements())) {
        m_tmpMatrix.reinit(m_systemMatrix->get_sparsity_pattern());
    }
#endif
    // A*f(t)
    m_tmpMatrix.copy_from(*m_systemMatrix);
    size_t n_blocks = m_systemMatrix->n_block_rows();
    cout << n_blocks << " Point 5" << endl;
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
;


//    dealii::SolverControl solver_control(1000, 1e-6*f.l2_norm(), false, false);	//* m_tmpSystemVector.l2_norm());
 //   dealii::SolverBicgstab<distributed_block_vector> bicgstab(solver_control);
#ifdef WITH_TRILINOS
  //  bicgstab.solve(m_tmpMatrix, result, f, dealii::PreconditionIdentity());


   // natrium::DealIIWrapper<distributed_sparse_block_matrix,
    // distributed_block_vector>::m_tmpMatrix = m_tmpMatrix;
	cout << m_tmpMatrix.memory_consumption() << endl;
    result = natrium::DealIIWrapper<distributed_sparse_block_matrix,distributed_block_vector>::
            solvers<dealii::PreconditionIdentity>(f, 1);
#else
    bicgstab.solve(m_tmpMatrix, result, f,
            dealii::PreconditionIdentity());	//,	           preconditioner);
#endif

    return result;

}

template<class MATRIX, class VECTOR>
double natrium::DealIIWrapper<MATRIX, VECTOR>::step(VECTOR& vector,
        const MATRIX& systemMatrix, const VECTOR& systemVector, double t,
        double dt) {

    if ((0.0 != dt) and dt != this->getTimeStepSize()){
        this->setTimeStepSize(dt);
        if (m_dealIIRKEmbedded){
            LOG(BASIC) << "Time step size set to " << dt << endl;
        } else {
            LOG(DETAILED) << "Time step size set to " << dt << endl;
        }
    } else if (0.0 == dt){
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
}

template<class MATRIX, class VECTOR>
natrium::DealIIWrapper<MATRIX, VECTOR>::DealIIWrapper(const double timeStepSize,
        const DealIntegratorName rkScheme, const DealSolverName linearSolver) :
        TimeIntegrator<MATRIX, VECTOR>(timeStepSize), m_systemMatrix(NULL), m_systemVector(
        NULL) {

    m_solver = linearSolver;
    assert(rkScheme < 12);
    assert(NONE != rkScheme);
    dealii::TimeStepping::runge_kutta_method rk =
            static_cast<dealii::TimeStepping::runge_kutta_method>((int) rkScheme);
    if (rkScheme < 3) {
        m_dealIIRKStepper = make_shared<
                dealii::TimeStepping::ExplicitRungeKutta<VECTOR> >(rk);
    } else if (rkScheme < 7) {
        m_dealIIRKStepper = make_shared<
                dealii::TimeStepping::ImplicitRungeKutta<VECTOR> >(rk);
    } else if (rkScheme < 12) {
        m_dealIIRKEmbedded = make_shared<
                dealii::TimeStepping::EmbeddedExplicitRungeKutta<VECTOR> >(rk);
        m_dealIIRKStepper = m_dealIIRKEmbedded;
    };
}

template<class MATRIX, class VECTOR>
natrium::DealIIWrapper<MATRIX, VECTOR>::DealIIWrapper(const double timeStepSize,
        const DealIntegratorName rkScheme, const DealSolverName linearSolver, double coarsen_param, double refine_param, double min_delta, double max_delta, double refine_tol, double coarsen_tol) :
        TimeIntegrator<MATRIX, VECTOR>(timeStepSize), m_systemMatrix(NULL), m_systemVector(
        NULL) {

    m_solver = linearSolver;


    assert(rkScheme < 12);
    assert(NONE != rkScheme);
    dealii::TimeStepping::runge_kutta_method rk =
            static_cast<dealii::TimeStepping::runge_kutta_method>((int) rkScheme);
    m_dealIIRKEmbedded = make_shared<
            dealii::TimeStepping::EmbeddedExplicitRungeKutta<VECTOR> >(rk,coarsen_param,refine_param,min_delta,max_delta,refine_tol,coarsen_tol);
    m_dealIIRKStepper = m_dealIIRKEmbedded;

}

/// explicit instantiation

template distributed_vector
        natrium::DealIIWrapper<distributed_sparse_matrix,distributed_vector>::solvers<distributed_sparse_matrix,distributed_vector,dealii::TrilinosWrappers::PreconditionIdentity>(const distributed_vector& f, bool isTmpMatrix);

template distributed_block_vector natrium::DealIIWrapper<distributed_sparse_block_matrix,distributed_block_vector>::
        solvers<distributed_sparse_block_matrix,distributed_block_vector,TrilinosBlockPreconditioner>(const distributed_block_vector& f, bool isTmpMatrix);

template distributed_block_vector natrium::DealIIWrapper<distributed_sparse_block_matrix,distributed_block_vector>::
        solvers<distributed_sparse_block_matrix,distributed_block_vector,dealii::PreconditionIdentity>(const distributed_block_vector& f, bool isTmpMatrix);

//template distributed_block_vector natrium::DealIIWrapper<distributed_sparse_block_matrix,distributed_block_vector>::
//        solve<dealii::TrilinosWrappers::PreconditionIdentity>(const distributed_block_vector& f) const;


template class DealIIWrapper<distributed_sparse_matrix, distributed_vector> ;
template class DealIIWrapper<distributed_sparse_block_matrix,
        distributed_block_vector> ;

} /* namespace natrium */
