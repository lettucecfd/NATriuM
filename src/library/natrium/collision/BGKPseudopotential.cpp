/*
 * BGKPseudopotential.cpp
 *
 *  Created on: Nov 5, 2014
 *      Author: kk
 */

#include "BGKPseudopotential.h"

#include <cmath>
#include "deal.II/dofs/dof_handler.h"

#include "deal.II/fe/fe_update_flags.h"
#include "../advection/SEDGMinLee.h"

#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"

#include <cassert>

#include "../utilities/Math.h"

namespace natrium {

/// constructor
BGKPseudopotential::BGKPseudopotential(double relaxationParameter, double dt,
		const shared_ptr<Stencil> stencil) :
		BGK(relaxationParameter, dt, stencil) {

	throw CollisionException("Pseudopotential model not implemented, yet.");

}

// destructor
BGKPseudopotential::~BGKPseudopotential() {
}

// getEquilibriumDistribution
double BGKPseudopotential::getEquilibriumDistribution(size_t i,
		const numeric_vector& u, const double rho) const {

	cout << "Function not implemented, yet." << endl;
	assert(false);

	assert(i < getStencil()->getQ());
	assert(rho > 0);
	assert(u.size() == getStencil()->getD());
	assert(u(0) < 1000000000000000.);
	assert(u(1) < 1000000000000000.);

	double prefactor = getStencil()->getWeight(i) * rho;
	double uSquareTerm = -Math::scalar_product(u, u)
			/ (2 * getStencil()->getSpeedOfSoundSquare());
	if (0 == i) {
		return prefactor * (1 + uSquareTerm);
	}
	double mixedTerm = Math::scalar_product(u, getStencil()->getDirection(i))
			/ getStencil()->getSpeedOfSoundSquare();
	return prefactor * (1 + mixedTerm * (1 + 0.5 * (mixedTerm)) + uSquareTerm);
}

// getInteractionForce
void BGKPseudopotential::getInteractionForce(
		const vector<double>& , numeric_vector & interactionForce,
		const double rho) {

	const double G = -4.4;
	numeric_vector densityGradient(2);

	interactionForce(0) = 0.0;
	interactionForce(1) = 0.0;

	//TODO getDensityGradient ;

	interactionForce(0) = -G * (1. - exp(rho)) * densityGradient(0);
	interactionForce(1) = -G * (1. - exp(rho)) * densityGradient(1);
}

/////////////////////////////
// Collision for PP model: //
/////////////////////////////
void BGKPseudopotential::collideAll(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		bool inInitializationProcedure,
		const dealii::IndexSet& locally_owned_dofs) const {

	if (Stencil_D2Q9 != getStencil()->getStencilType()) {
		// Inefficient collision for other than D2Q9
		BGK::collideAll(f, densities, velocities, inInitializationProcedure);
	} else {
		size_t n_dofs = f.at(0).size();
		size_t Q = 9;
		size_t D = 2;
		double scaling = getStencil()->getScaling();
		double cs2 = getStencil()->getSpeedOfSoundSquare();
		double prefactor = scaling / cs2;
		double relax_factor = getPrefactor();

		assert(f.size() == Q);
		assert(velocities.size() == D);

#ifdef DEBUG
		for (size_t i = 0; i < Q; i++) {
			assert (f.at(i).size() == n_dofs);
		}
		assert (densities.size() == n_dofs);
		for (size_t i = 0; i < D; i++) {
			assert (velocities.at(i).size() == n_dofs);
		}
#endif

		// allocation
		vector<double> feq(Q, 0.0);
		double scalar_product;
		double uSquareTerm;
		double mixedTerm;
		double weighting;

		// for all cells
		// TODO (when fully implemented) profiling! is the cost of this procedure OK (especially concerning cache exploitation)
		// TODO (KNUT) The PP_Model class needs a (shared) pointer to AdvectionOperator and a function getAdvectionOperator()
		const dealii::DoFHandler<2>& dof_handler =
				*m_advectionOperator->getDoFHandler();
		const size_t dofs_per_cell =
				m_advectionOperator->getNumberOfDoFsPerCell();
		std::vector<dealii::types::global_dof_index> localDoFIndices(
				dofs_per_cell);
		const std::map<size_t, size_t>& celldofToQIndex =
				m_advectionOperator->getCelldofToQIndex();
		const dealii::UpdateFlags updateFlags = dealii::update_gradients;

		dealii::FEValues<2> feValues(m_advectionOperator->getMapping(),
				*m_advectionOperator->getFe(),
				*m_advectionOperator->getQuadrature(), updateFlags);

		// loop over all cells
		dealii::Tensor<1, 2> density_gradient;
		typename dealii::DoFHandler<2>::active_cell_iterator cell =
				dof_handler.begin_active(), endc = dof_handler.end();
		for (; cell != endc; ++cell) {

			if (! cell->is_locally_owned())
				continue;

			// get global degrees of freedom
			cell->get_dof_indices(localDoFIndices);

			// calculate gradients of shape functions
			feValues.reinit(cell);

			// for all dofs in cell
			for (size_t j = 0; j < dofs_per_cell; j++) {

				size_t i = localDoFIndices.at(j);

				// calculate density
				densities(i) = f.at(0)(i) + f.at(1)(i) + f.at(2)(i) + f.at(3)(i)
						+ f.at(4)(i) + f.at(5)(i) + f.at(6)(i) + f.at(7)(i)
						+ f.at(8)(i);
				if (densities(i) < 1e-10) {
					throw CollisionException(
							"Densities too small (< 1e-10) for collisions. Decrease time step size.");
				}

				// calculate density gradient
				// This function assumes that the DoFs are uniquely identified with quadrature points
				// it uses the relation grad(rho)(x_i) = grad( sum_alpha f_alpha(x_i))  = grad ( sum_alpha(sum_k(dof_{alpha,k} phi_k(x_i))))
				//                                     = sum_k ((sum_alpha dof_{alpha,k}) grad(phi_k(x_i))) = sum_k (rho_k grad(phi_k(x_i)))
				for (size_t k = 0; k < dofs_per_cell; k++) {
					// TODO: make this line work
				//	density_gradient += densities(localDoFIndices.at(k))
				//			* feValues.shape_grad(k, celldofToQIndex.at(j));
				}

				// iterative initialization scheme (solve poisson equation indirectly)
				if (not inInitializationProcedure) {
					// calculate velocity
					// for all velocity components
					velocities.at(0)(i) = scaling / densities(i)
							* (f.at(1)(i) + f.at(5)(i) + f.at(8)(i) - f.at(3)(i)
									- f.at(6)(i) - f.at(7)(i));
					velocities.at(1)(i) = scaling / densities(i)
							* (f.at(2)(i) + f.at(5)(i) + f.at(6)(i) - f.at(4)(i)
									- f.at(7)(i) - f.at(8)(i));
				}

				// TODO (KNUT) Change to PP Model!
				// calculate equilibrium distribution
				scalar_product = velocities.at(0)(i) * velocities.at(0)(i)
						+ velocities.at(1)(i) * velocities.at(1)(i);
				uSquareTerm = -scalar_product / (2 * cs2);
				// direction 0
				weighting = 4. / 9. * densities(i);
				feq.at(0) = weighting * (1 + uSquareTerm);
				// directions 1-4
				weighting = 1. / 9. * densities(i);
				mixedTerm = prefactor * (velocities.at(0)(i));
				feq.at(1) = weighting
						* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
				feq.at(3) = weighting
						* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
				mixedTerm = prefactor * (velocities.at(1)(i));
				feq.at(2) = weighting
						* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
				feq.at(4) = weighting
						* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
				// directions 5-8
				weighting = 1. / 36. * densities(i);
				mixedTerm = prefactor
						* (velocities.at(0)(i) + velocities.at(1)(i));
				feq.at(5) = weighting
						* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
				feq.at(7) = weighting
						* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
				mixedTerm = prefactor
						* (-velocities.at(0)(i) + velocities.at(1)(i));
				feq.at(6) = weighting
						* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
				feq.at(8) = weighting
						* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);

				// BGK collision
				f.at(0)(i) += relax_factor * (f.at(0)(i) - feq.at(0));
				f.at(1)(i) += relax_factor * (f.at(1)(i) - feq.at(1));
				f.at(2)(i) += relax_factor * (f.at(2)(i) - feq.at(2));
				f.at(3)(i) += relax_factor * (f.at(3)(i) - feq.at(3));
				f.at(4)(i) += relax_factor * (f.at(4)(i) - feq.at(4));
				f.at(5)(i) += relax_factor * (f.at(5)(i) - feq.at(5));
				f.at(6)(i) += relax_factor * (f.at(6)(i) - feq.at(6));
				f.at(7)(i) += relax_factor * (f.at(7)(i) - feq.at(7));
				f.at(8)(i) += relax_factor * (f.at(8)(i) - feq.at(8));
			} /* for all dofs in cell */

		} /* for all cells */

	} /* if D2Q9 */
}

}
