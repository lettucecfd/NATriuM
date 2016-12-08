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
double n[][9] = { { 8, 2, 2, 2, 2, -4, -4, -4, -4 }, { 1, 13, -1, 1,
		-1, 5, -1, -1, 5 }, { 1, -1, 13, -1, 1, 5, 5, -1, -1 }, { 1, 1, -1,
		13, -1, -1, 5, 5, -1 }, { 1, -1, 1, -1, 13, -1, -1, 5, 5 }, { -1, 5,
		5, -1, -1, 5, -1, 2, -1 }, { -1, -1, 5, 5, -1, -1, 5, -1, 2 }, { -1,
		-1, -1, 5, 5, 2, -1, 5, -1 }, { -1, 5, -1, -1, 5, -1, 2, -1, 5 } };
// denominator
double d[][9] = { { 9, 9, 9, 9, 9, 9, 9, 9, 9 }, { 18, 18, 9, 18, 9,
		9, 9, 9, 9 }, { 18, 9, 18, 9, 18, 9, 9, 9, 9 }, { 18, 18, 9, 18, 9,
		9, 9, 9, 9 }, { 18, 9, 18, 9, 18, 9, 9, 9, 9 }, { 36, 36, 36, 36,
		36, 9, 9, 9, 9 }, { 36, 36, 36, 36, 36, 9, 9, 9, 9 }, { 36, 36, 36,
		36, 36, 9, 9, 9, 9 }, { 36, 36, 36, 36, 36, 9, 9, 9, 9 } };


template<size_t dim>
void PseudoEntropicStabilizer<dim>::apply() {


	TimerOutput::Scope timer_section(Timing::getTimer(), "Pseudo-entropic stabilizer");

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

		// calculate conserved
		double rho = 0;
		double jx = 0;
		double jy = 0;
		for (size_t iQ = 0; iQ < stencil.getQ(); ++iQ) {
			rho += f_i[iQ];
			jx += f_i[iQ] * stencil.getDirection(iQ)[0];
			jy += f_i[iQ] * stencil.getDirection(iQ)[1];
		}


		for (size_t i = 0; i < 9; i++) {
			f_i_new[i] = 0;
			for (size_t j = 0; j < 9; j++) {
				f_i_new[i] += n[i][j] / d[i][j] * f_i[j];
			}
		}

		double rho_pc = 0;
		double jx_pc = 0;
		double jy_pc = 0;
		for (size_t iQ = 0; iQ < stencil.getQ(); ++iQ) {
			rho_pc += f_i_new[iQ];
			jx_pc += f_i_new[iQ] * stencil.getDirection(iQ)[0];
			jy_pc += f_i_new[iQ] * stencil.getDirection(iQ)[1];
		}

		assert (fabs(rho_pc - rho) < 1e-8);
		assert (fabs(jx_pc - jx) < 1e-8);
		assert (fabs(jy_pc - jy) < 1e-8);

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
} /* apply */

template class PseudoEntropicStabilizer<2> ;
template class PseudoEntropicStabilizer<3> ;

} /* namespace natrium */
