/**
 * @file CollisionModel.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "Collision.h"

#include <cmath>

#include "deal.II/dofs/dof_handler.h"

#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"
#include "deal.II/fe/fe_update_flags.h"

#include "../solver/DistributionFunctions.h"

#include "../boltzmannmodels/D2Q9IncompressibleModel.h"
#include "../boltzmannmodels/D2Q9PseudopotentialModel.h"

#include "../advection/SEDGMinLee.h"

#include "BGKTransformed.h"

namespace natrium {

// constructor
template<class BoltzmannType, class CollisionType>
Collision<BoltzmannType, CollisionType>::Collision(BoltzmannType boltzmannType,
		CollisionType collisionType) :
		m_boltzmannType(boltzmannType), m_collisionType(collisionType), m_d(
				boltzmannType.getD()), m_q(boltzmannType.getQ()) {

} // constructor
template Collision<D2Q9IncompressibleModel, BGKTransformed>::Collision(
		D2Q9IncompressibleModel boltzmannType, BGKTransformed collisionModel);
template Collision<D2Q9PseudopotentialModel, BGKTransformed>::Collision(
		D2Q9PseudopotentialModel boltzmannType, BGKTransformed collisionModel);


template<class BoltzmannType, class CollisionType>
Collision<BoltzmannType, CollisionType>::~Collision() {
}
template Collision<D2Q9IncompressibleModel, BGKTransformed>::~Collision();
template Collision<D2Q9PseudopotentialModel, BGKTransformed>::~Collision();

template<class BoltzmannType, class CollisionType>
void Collision<BoltzmannType, CollisionType>::collideAll(
		DistributionFunctions& f, distributed_vector& densities,
		vector<distributed_vector>& velocities,
		bool inInitializationProcedure) const {

	size_t n_dofs = f.at(0).size();
	size_t Q = m_boltzmannType.getQ();

	assert(f.size() == Q);
	assert(velocities.size() == m_d);

#ifdef DEBUG
	for (size_t i = 0; i < Q; i++) {
		assert (f.at(i).size() == n_dofs);
	}
	assert (densities.size() == n_dofs);
	for (size_t i = 0; i < m_d; i++) {
		assert (velocities.at(i).size() == n_dofs);
	}
#endif

//for all degrees of freedom
	for (size_t i = 0; i < n_dofs; i++) {

		// calculate density
		densities(i) = 0;
		for (size_t j = 0; j < Q; j++) {
			densities(i) += f.at(j)(i);
		}
		if (densities(i) < 1e-10) {
			throw CollisionException(
					"Densities too small (< 1e-10) for collisions. Decrease time step size.");
		}

		if (not inInitializationProcedure) {
			// calculate velocity
			// for all velocity components
			for (size_t j = 0; j < m_d; j++) {
				velocities.at(j)(i) = 0;

				for (size_t k = 0; k < Q; k++) {
					velocities.at(j)(i) += f.at(k)(i)
							* m_boltzmannType.getDirection(k)(j);
				}
				velocities.at(j)(i) /= densities(i);
			}
		}

		// calculate equilibrium distribution
		// TODO Optimize by passing different arguments
		vector<double> feq(Q, 0.0);
		numeric_vector u(m_d);
		for (size_t j = 0; j < m_d; j++) {
			u(j) = velocities.at(j)(i);
		}
		m_boltzmannType.getEquilibriumDistributions(feq, u, densities(i));

		// BGK collision
		m_collisionType.collideSingleDoF(i, feq, f);
	}
}
template<> void Collision<D2Q9IncompressibleModel, BGKTransformed>::collideAll(
		DistributionFunctions& f, distributed_vector& densities,
		vector<distributed_vector>& velocities,
		bool inInitializationProcedure) const {

	size_t n_dofs = f.at(0).size();
	size_t Q = 9;
	size_t D = 2;
	double scaling = m_boltzmannType.getScaling();
	double cs2 = m_boltzmannType.getSpeedOfSoundSquare();
	double prefactor = scaling / cs2;
	double relax_factor = m_collisionType.getPrefactor();

	assert(f.size() == Q);
	assert(velocities.size() == D);

#ifdef DEBUG
	for (size_t i = 0; i < Q; i++) {
		assert (f.at(i).size() == n_dofs);
	}
	assert (densities.size() == n_dofs);
	for (size_t i = 0; i < m_d; i++) {
		assert (velocities.at(i).size() == n_dofs);
	}
#endif

	// allocation
	vector<double> feq(Q, 0.0);
	double scalar_product;
	double uSquareTerm;
	double mixedTerm;
	double weighting;

	// for all dofs
	for (size_t i = 0; i < n_dofs; i++) {

		// calculate density
		densities(i) = f.at(0)(i) + f.at(1)(i) + f.at(2)(i) + f.at(3)(i)
				+ f.at(4)(i) + f.at(5)(i) + f.at(6)(i) + f.at(7)(i)
				+ f.at(8)(i);
		if (densities(i) < 1e-10) {
			throw CollisionException(
					"Densities too small (< 1e-10) for collisions. Decrease time step size.");
		}

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
		mixedTerm = prefactor * (velocities.at(0)(i) + velocities.at(1)(i));
		feq.at(5) = weighting
				* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
		feq.at(7) = weighting
				* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
		mixedTerm = prefactor * (-velocities.at(0)(i) + velocities.at(1)(i));
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
	}
}

/////////////////////////////
// Collision for PP model: //
/////////////////////////////
template<> void Collision<D2Q9PseudopotentialModel, BGKTransformed>::collideAll(
		DistributionFunctions& f, distributed_vector& densities,
		vector<distributed_vector>& velocities,
		bool inInitializationProcedure) const {
	size_t n_dofs = f.at(0).size();
	size_t Q = 9;
	size_t D = 2;
	double scaling = m_boltzmannType.getScaling();
	double cs2 = m_boltzmannType.getSpeedOfSoundSquare();
	double prefactor = scaling / cs2;
	double relax_factor = m_collisionType.getPrefactor();

	assert(f.size() == Q);
	assert(velocities.size() == D);

#ifdef DEBUG
	for (size_t i = 0; i < Q; i++) {
		assert (f.at(i).size() == n_dofs);
	}
	assert (densities.size() == n_dofs);
	for (size_t i = 0; i < m_d; i++) {
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
	const SEDGMinLee<2>& advectionOperator = *m_boltzmannType.getAdvectionOperator();
	const dealii::DoFHandler<2>& dof_handler = *advectionOperator.getDoFHandler();
	const size_t dofs_per_cell = advectionOperator.getNumberOfDoFsPerCell();
	std::vector<dealii::types::global_dof_index> localDoFIndices(dofs_per_cell);
	const std::map<size_t, size_t>& celldofToQIndex = advectionOperator.getCelldofToQIndex();
	const dealii::UpdateFlags updateFlags = dealii::update_gradients;

	dealii::FEValues<2> feValues(advectionOperator.getMapping(), *advectionOperator.getFe(), *advectionOperator.getQuadrature(), updateFlags);



	// loop over all cells
	dealii::Tensor<1,2> density_gradient;
	typename dealii::DoFHandler<2>::active_cell_iterator cell =
			dof_handler.begin_active(), endc = dof_handler.end();
	for (; cell != endc; ++cell) {

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
			for (size_t k = 0; k < dofs_per_cell; k++){
				density_gradient += densities(localDoFIndices.at(k))*feValues.shape_grad(k, celldofToQIndex.at(j));
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
			mixedTerm = prefactor * (velocities.at(0)(i) + velocities.at(1)(i));
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

}

template<class BoltzmannType, class CollisionType>
void Collision<BoltzmannType, CollisionType>::collideAll(
		DistributionFunctions& f, distributed_vector& densities,
		vector<distributed_vector>& velocities, const vector<bool>& isBoundary,
		bool inInitializationProcedure) const {

	size_t n_dofs = f.at(0).size();
	size_t Q = m_boltzmannType.getQ();

	assert(f.size() == Q);
	assert(velocities.size() == m_d);

#ifdef DEBUG
	for (size_t i = 0; i < Q; i++) {
		assert (f.at(i).size() == n_dofs);
	}
	assert (densities.size() == n_dofs);
	for (size_t i = 0; i < m_d; i++) {
		assert (velocities.at(i).size() == n_dofs);
	}
#endif

//for all degrees of freedom
	for (size_t i = 0; i < n_dofs; i++) {

		// calculate density
		densities(i) = 0;
		for (size_t j = 0; j < Q; j++) {
			densities(i) += f.at(j)(i);
		}
		if (densities(i) < 1e-10) {
			throw CollisionException(
					"Densities too small (< 1e-10) for collisions. Decrease time step size.");
		}

		if (not inInitializationProcedure) {
			// calculate velocity
			// for all velocity components
			for (size_t j = 0; j < m_d; j++) {
				velocities.at(j)(i) = 0;

				for (size_t k = 0; k < Q; k++) {
					velocities.at(j)(i) += f.at(k)(i)
							* m_boltzmannType.getDirection(k)(j);
				}
				velocities.at(j)(i) /= densities(i);
			}
		}

		if (isBoundary.at(i)) {
			continue;
		}

		// calculate equilibrium distribution
		// TODO Optimize by passing different arguments
		vector<double> feq(Q, 0.0);
		numeric_vector u(m_d);
		for (size_t j = 0; j < m_d; j++) {
			u(j) = velocities.at(j)(i);
		}
		m_boltzmannType.getEquilibriumDistributions(feq, u, densities(i));

		// BGK collision
		m_collisionType.collideSingleDoF(i, feq, f);
	}
} /*collideAll */
template void Collision<D2Q9IncompressibleModel, BGKTransformed>::collideAll(
		DistributionFunctions& f, distributed_vector& densities,
		vector<distributed_vector>& velocities, const vector<bool>& isBoundary,
		bool inInitializationProcedure) const;

} /* namespace natrium */
