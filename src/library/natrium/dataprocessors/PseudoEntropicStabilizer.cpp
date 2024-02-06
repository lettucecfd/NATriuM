/**
 * @file BGKStandard.cpp
 * @short D2Q9 model description for incompressible flow.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "PseudoEntropicStabilizer.h"
#include "../stencils/Stencil.h"
#include "../solver/CFDSolver.h"
#include "../solver/DistributionFunctions.h"
#include "deal.II/dofs/dof_handler.h"
#include "deal.II/base/index_set.h"
#include "../utilities/Timing.h"

// enable + operator for filling vectors
#include "boost/assign/std/vector.hpp"

// enable + operator for filling vectors
using namespace boost::assign;

namespace natrium {

//dealii::FullMatrix stabilizer_matrix(9);
// nominator
double n[][9] = { { 8, 2, 2, 2, 2, -4, -4, -4, -4 }, { 1, 13, -1, 1, -1, 5, -1,
		-1, 5 }, { 1, -1, 13, -1, 1, 5, 5, -1, -1 }, { 1, 1, -1, 13, -1, -1, 5,
		5, -1 }, { 1, -1, 1, -1, 13, -1, -1, 5, 5 }, { -1, 5, 5, -1, -1, 5, -1,
		2, -1 }, { -1, -1, 5, 5, -1, -1, 5, -1, 2 }, { -1, -1, -1, 5, 5, 2, -1,
		5, -1 }, { -1, 5, -1, -1, 5, -1, 2, -1, 5 } };
// denominator
double d[][9] = { { 9, 9, 9, 9, 9, 9, 9, 9, 9 },
		{ 18, 18, 9, 18, 9, 9, 9, 9, 9 }, { 18, 9, 18, 9, 18, 9, 9, 9, 9 }, {
				18, 18, 9, 18, 9, 9, 9, 9, 9 },
		{ 18, 9, 18, 9, 18, 9, 9, 9, 9 }, { 36, 36, 36, 36, 36, 9, 9, 9, 9 }, {
				36, 36, 36, 36, 36, 9, 9, 9, 9 }, { 36, 36, 36, 36, 36, 9, 9, 9,
				9 }, { 36, 36, 36, 36, 36, 9, 9, 9, 9 } };

double nd_d2q9_with_e[][9] = { { 4. / 9., 4. / 9., 4. / 9., 4. / 9., 4. / 9., 4.
		/ 9., 4. / 9., 4. / 9., 4. / 9. }, { 1. / 9., 25. / 36., -5. / 36., 1.
		/ 36., -5. / 36., 4. / 9., -2. / 9., -2. / 9., 4. / 9. }, { 1. / 9., -5.
		/ 36., 25. / 36., -5. / 36., 1. / 36., 4. / 9., 4. / 9., -2. / 9., -2.
		/ 9. }, { 1. / 9., 1. / 36., -5. / 36., 25. / 36., -5. / 36., -2. / 9.,
		4. / 9., 4. / 9., -2. / 9. }, { 1. / 9., -5. / 36., 1. / 36., -5. / 36.,
		25. / 36., -2. / 9., -2. / 9., 4. / 9., 4. / 9. }, { 1. / 36., 1. / 9.,
		1. / 9., -1. / 18., -1. / 18., 4. / 9., -2. / 9., 1. / 9., -2. / 9. }, {
		1. / 36., -1. / 18., 1. / 9., 1. / 9., -1. / 18., -2. / 9., 4. / 9., -2.
				/ 9., 1. / 9. }, { 1. / 36., -1. / 18., -1. / 18., 1. / 9., 1.
		/ 9., 1. / 9., -2. / 9., 4. / 9., -2. / 9. }, { 1. / 36., 1. / 9., -1.
		/ 18., -1. / 18., 1. / 9., -2. / 9., 1. / 9., -2. / 9., 4. / 9. } };

double nd_d3q19[][19] = { { 5. / 6., 1. / 3., 1. / 3., 1. / 3., 1. / 3., 1.
		/ 3., 1. / 3., -1. / 6., -1. / 6., -1. / 6., -1. / 6., -1. / 6., -1.
		/ 6., -1. / 6., -1. / 6., -1. / 6., -1. / 6., -1. / 6., -1. / 6. }, { 1.
		/ 18., 7. / 18., -1. / 36., 1. / 18., -1. / 36., -1. / 36., -1. / 36.,
		11. / 36., -1. / 36., -1. / 36., 11. / 36., 11. / 36., 11. / 36., -1.
				/ 36., -1. / 36., -1. / 9., -1. / 9., -1. / 9., -1. / 9. }, { 1.
		/ 18., -1. / 36., 7. / 18., -1. / 36., 1. / 18., -1. / 36., -1. / 36.,
		11. / 36., 11. / 36., -1. / 36., -1. / 36., -1. / 9., -1. / 9., -1.
				/ 9., -1. / 9., 11. / 36., 11. / 36., -1. / 36., -1. / 36. },
		{ 1. / 18., 1. / 18., -1. / 36., 7. / 18., -1. / 36., -1. / 36., -1.
				/ 36., -1. / 36., 11. / 36., 11. / 36., -1. / 36., -1. / 36.,
				-1. / 36., 11. / 36., 11. / 36., -1. / 9., -1. / 9., -1. / 9.,
				-1. / 9. }, { 1. / 18., -1. / 36., 1. / 18., -1. / 36., 7.
				/ 18., -1. / 36., -1. / 36., -1. / 36., -1. / 36., 11. / 36.,
				11. / 36., -1. / 9., -1. / 9., -1. / 9., -1. / 9., -1. / 36.,
				-1. / 36., 11. / 36., 11. / 36. }, { 1. / 18., -1. / 36., -1.
				/ 36., -1. / 36., -1. / 36., 7. / 18., 1. / 18., -1. / 9., -1.
				/ 9., -1. / 9., -1. / 9., 11. / 36., -1. / 36., -1. / 36., 11.
				/ 36., 11. / 36., -1. / 36., -1. / 36., 11. / 36. },
		{ 1. / 18., -1. / 36., -1. / 36., -1. / 36., -1. / 36., 1. / 18., 7.
				/ 18., -1. / 9., -1. / 9., -1. / 9., -1. / 9., -1. / 36., 11.
				/ 36., 11. / 36., -1. / 36., -1. / 36., 11. / 36., 11. / 36.,
				-1. / 36. }, { -1. / 72., 11. / 72., 11. / 72., -1. / 72., -1.
				/ 72., -1. / 18., -1. / 18., 41. / 72., -7. / 72., 17. / 72.,
				-7. / 72., 1. / 9., 1. / 9., -1. / 18., -1. / 18., 1. / 9., 1.
						/ 9., -1. / 18., -1. / 18. }, { -1. / 72., -1. / 72.,
				11. / 72., 11. / 72., -1. / 72., -1. / 18., -1. / 18., -7.
						/ 72., 41. / 72., -7. / 72., 17. / 72., -1. / 18., -1.
						/ 18., 1. / 9., 1. / 9., 1. / 9., 1. / 9., -1. / 18.,
				-1. / 18. }, { -1. / 72., -1. / 72., -1. / 72., 11. / 72., 11.
				/ 72., -1. / 18., -1. / 18., 17. / 72., -7. / 72., 41. / 72.,
				-7. / 72., -1. / 18., -1. / 18., 1. / 9., 1. / 9., -1. / 18.,
				-1. / 18., 1. / 9., 1. / 9. }, { -1. / 72., 11. / 72., -1.
				/ 72., -1. / 72., 11. / 72., -1. / 18., -1. / 18., -7. / 72.,
				17. / 72., -7. / 72., 41. / 72., 1. / 9., 1. / 9., -1. / 18.,
				-1. / 18., -1. / 18., -1. / 18., 1. / 9., 1. / 9. }, { -1.
				/ 72., 11. / 72., -1. / 18., -1. / 72., -1. / 18., 11. / 72.,
				-1. / 72., 1. / 9., -1. / 18., -1. / 18., 1. / 9., 41. / 72.,
				-7. / 72., 17. / 72., -7. / 72., 1. / 9., -1. / 18., -1. / 18.,
				1. / 9. }, { -1. / 72., 11. / 72., -1. / 18., -1. / 72., -1.
				/ 18., -1. / 72., 11. / 72., 1. / 9., -1. / 18., -1. / 18., 1.
				/ 9., -7. / 72., 41. / 72., -7. / 72., 17. / 72., -1. / 18., 1.
				/ 9., 1. / 9., -1. / 18. }, { -1. / 72., -1. / 72., -1. / 18.,
				11. / 72., -1. / 18., -1. / 72., 11. / 72., -1. / 18., 1. / 9.,
				1. / 9., -1. / 18., 17. / 72., -7. / 72., 41. / 72., -7. / 72.,
				-1. / 18., 1. / 9., 1. / 9., -1. / 18. }, { -1. / 72., -1.
				/ 72., -1. / 18., 11. / 72., -1. / 18., 11. / 72., -1. / 72.,
				-1. / 18., 1. / 9., 1. / 9., -1. / 18., -7. / 72., 17. / 72.,
				-7. / 72., 41. / 72., 1. / 9., -1. / 18., -1. / 18., 1. / 9. },
		{ -1. / 72., -1. / 18., 11. / 72., -1. / 18., -1. / 72., 11. / 72., -1.
				/ 72., 1. / 9., 1. / 9., -1. / 18., -1. / 18., 1. / 9., -1.
				/ 18., -1. / 18., 1. / 9., 41. / 72., -7. / 72., 17. / 72., -7.
				/ 72. }, { -1. / 72., -1. / 18., 11. / 72., -1. / 18., -1.
				/ 72., -1. / 72., 11. / 72., 1. / 9., 1. / 9., -1. / 18., -1.
				/ 18., -1. / 18., 1. / 9., 1. / 9., -1. / 18., -7. / 72., 41.
				/ 72., -7. / 72., 17. / 72. }, { -1. / 72., -1. / 18., -1.
				/ 72., -1. / 18., 11. / 72., -1. / 72., 11. / 72., -1. / 18.,
				-1. / 18., 1. / 9., 1. / 9., -1. / 18., 1. / 9., 1. / 9., -1.
						/ 18., 17. / 72., -7. / 72., 41. / 72., -7. / 72. }, {
				-1. / 72., -1. / 18., -1. / 72., -1. / 18., 11. / 72., 11.
						/ 72., -1. / 72., -1. / 18., -1. / 18., 1. / 9., 1.
						/ 9., 1. / 9., -1. / 18., -1. / 18., 1. / 9., -7. / 72.,
				17. / 72., -7. / 72., 41. / 72. } };

template<size_t dim>
void PseudoEntropicStabilizer<dim>::apply_d2q9() {

	//pout << "Stabilizer is active." << endl;
	const Stencil& stencil = *(this->m_solver.getStencil());
	assert(Stencil_D2Q9 == stencil.getStencilType());

	// allocation
	double f_i[9] = { };
	double f_i_new[9] = { };

	const dealii::IndexSet& locally_owned_dofs =
			this->m_solver.getAdvectionOperator()->getDoFHandler()->locally_owned_dofs();
	DistributionFunctions& f = this->m_solver.m_f;

	//for all degrees of freedom on current processor
	dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
	dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());
	for (it = locally_owned_dofs.begin(); it != end; it++) {
		size_t i = *it;

		// copy fi
		f_i[0] = f.at(0)(i);
		f_i[1] = f.at(1)(i);
		f_i[2] = f.at(2)(i);
		f_i[3] = f.at(3)(i);
		f_i[4] = f.at(4)(i);
		f_i[5] = f.at(5)(i);
		f_i[6] = f.at(6)(i);
		f_i[7] = f.at(7)(i);
		f_i[8] = f.at(8)(i);

		for (size_t i = 0; i < 9; i++) {
			f_i_new[i] = 0;
			for (size_t j = 0; j < 9; j++) {
				f_i_new[i] += n[i][j] / d[i][j] * f_i[j];
			}
		}

		// copy to global variables
		f.at(0)(i) = f_i_new[0];
		f.at(1)(i) = f_i_new[1];
		f.at(2)(i) = f_i_new[2];
		f.at(3)(i) = f_i_new[3];
		f.at(4)(i) = f_i_new[4];
		f.at(5)(i) = f_i_new[5];
		f.at(6)(i) = f_i_new[6];
		f.at(7)(i) = f_i_new[7];
		f.at(8)(i) = f_i_new[8];

	} /* for all dofs */

	this->m_solver.m_f.updateGhosted();
} /* apply d2q9 */

template<size_t dim>
void PseudoEntropicStabilizer<dim>::apply_d2q9_with_e() {

	//pout << "Stabilizer is active." << endl;
	const Stencil& stencil = *(this->m_solver.getStencil());
	assert(Stencil_D2Q9 == stencil.getStencilType());

	// allocation
	double f_i[9] = { };
	double f_i_new[9] = { };

	const dealii::IndexSet& locally_owned_dofs =
			this->m_solver.getAdvectionOperator()->getDoFHandler()->locally_owned_dofs();
	DistributionFunctions& f = this->m_solver.m_f;

	//for all degrees of freedom on current processor
	dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
	dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());
	for (it = locally_owned_dofs.begin(); it != end; it++) {
		size_t i = *it;

		// copy fi
		f_i[0] = f.at(0)(i);
		f_i[1] = f.at(1)(i);
		f_i[2] = f.at(2)(i);
		f_i[3] = f.at(3)(i);
		f_i[4] = f.at(4)(i);
		f_i[5] = f.at(5)(i);
		f_i[6] = f.at(6)(i);
		f_i[7] = f.at(7)(i);
		f_i[8] = f.at(8)(i);

		for (size_t i = 0; i < 9; i++) {
			f_i_new[i] = 0;
			for (size_t j = 0; j < 9; j++) {
				f_i_new[i] += nd_d2q9_with_e[i][j] * f_i[j];
			}
		}

		// copy to global variables
		f.at(0)(i) = f_i_new[0];
		f.at(1)(i) = f_i_new[1];
		f.at(2)(i) = f_i_new[2];
		f.at(3)(i) = f_i_new[3];
		f.at(4)(i) = f_i_new[4];
		f.at(5)(i) = f_i_new[5];
		f.at(6)(i) = f_i_new[6];
		f.at(7)(i) = f_i_new[7];
		f.at(8)(i) = f_i_new[8];

	} /* for all dofs */

	this->m_solver.m_f.updateGhosted();
} /* apply d2q9_with_e */

template<size_t dim>
void PseudoEntropicStabilizer<dim>::apply_d3q19() {

	//pout << "Stabilizer is active." << endl;
	const Stencil& stencil = *(this->m_solver.getStencil());
	assert(Stencil_D3Q19 == stencil.getStencilType());

	// allocation
	double f_i[19] = { };
	double f_i_new[19] = { };

	const dealii::IndexSet& locally_owned_dofs =
			this->m_solver.getAdvectionOperator()->getDoFHandler()->locally_owned_dofs();
	DistributionFunctions& f = this->m_solver.m_f;

	//for all degrees of freedom on current processor
	dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
	dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());
	for (it = locally_owned_dofs.begin(); it != end; it++) {
		size_t i = *it;

		// copy fi
		for (size_t j = 0; j < 19; j++) {
			f_i[j] = f.at(j)(i);
		}

		// multiply stabilizer matrix
		for (size_t j = 0; j < 19; j++) {
			f_i_new[j] = 0;
			for (size_t k = 0; k < 19; k++) {
				f_i_new[j] += nd_d3q19[j][k] * f_i[k];
			}
		}

		// copy to global variables
		for (size_t j = 0; j < 19; j++) {
			f.at(j)(i) = f_i_new[j];
		}

	} /* for all dofs */

	this->m_solver.m_f.updateGhosted();
} /* apply d3q19 */

template<size_t dim>
void PseudoEntropicStabilizer<dim>::apply() {

	TimerOutput::Scope timer_section(Timing::getTimer(),
			"Pseudo-entropic stabilizer");

	//pout << "Stabilizer is active." << endl;
	const Stencil& stencil = *(this->m_solver.getStencil());

	if (Stencil_D2Q9 == stencil.getStencilType()) {
		if (m_withE)
			apply_d2q9_with_e();
		else
			apply_d2q9();
	} else if (Stencil_D3Q19 == stencil.getStencilType()) {
		apply_d3q19();
	}

	this->m_solver.m_f.updateGhosted();
} /* apply */




template<> void applyStabilizer<9>(const array<double,9>& in, array<double,9>& out){
//    assert( dim == 2);
	// multiply stabilizer matrix
	for (size_t j = 0; j < 9; j++) {
		out[j] = 0;
		for (size_t k = 0; k < 9; k++) {
			out[j] += n[j][k] / d[j][k] * in[k];
		}
	}
}

template<> void applyStabilizer<19>(const array<double,19>& in, array<double,19>& out){
    //assert( dim == 3);
	// multiply stabilizer matrix
	for (size_t j = 0; j < 19; j++) {
		out[j] = 0;
		for (size_t k = 0; k < 19; k++) {
			 out[j] += nd_d3q19[j][k] * in[k];
		}
	}
}

template<> void applyStabilizer<27>(const array<double,27>& in, array<double,27>& out){
//    assert( dim == 3);
    (void)in;
    (void)out;
    cout << "template<> void applyStabilizer<27> NOT IMPLEMENTED, YET." << endl;
    assert (false);
}

template class PseudoEntropicStabilizer<2> ;
template class PseudoEntropicStabilizer<3> ;

} /* namespace natrium */
