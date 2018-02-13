/*
 * BGKPseudopotential.cpp
 *
 *  Created on: Nov 5, 2014
 *      Author: kk
 */

#include "BGKPseudopotential.h"

#include <cmath>
#include <cassert>

#include "deal.II/dofs/dof_handler.h"
#include "deal.II/fe/fe_update_flags.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"

#include "ExternalForceFunctions.h"
#include "../utilities/Math.h"
#include "../advection/SEDGMinLee.h"

namespace natrium {

/// constructor
template<size_t dim>
BGKPseudopotential<dim>::BGKPseudopotential(double relaxationParameter,
		double dt, const boost::shared_ptr<Stencil> stencil,
		PseudopotentialParameters parameters) :
		BGK(relaxationParameter, dt, stencil), m_pseudopotentialParameters(
				parameters) {
	assert(parameters.G != 0);
	if (m_pseudopotentialParameters.pseudopotentialType == SHAN_CHEN) {

	} else if (m_pseudopotentialParameters.pseudopotentialType == SUKOP) {

	} else if (m_pseudopotentialParameters.pseudopotentialType
			== CARNAHAN_STARLING) {
		assert(parameters.T != 0);
	}

}

// destructor
template<size_t dim>
BGKPseudopotential<dim>::~BGKPseudopotential() {
}

// getEquilibriumDistribution
template<size_t dim>
double BGKPseudopotential<dim>::getEquilibriumDistribution(size_t i,
		const numeric_vector& u, const double rho) const {

	//pout << "Function not implemented, yet." << endl;
	//assert(false);

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
template<size_t dim>
void BGKPseudopotential<dim>::getInteractionForce(const vector<double>&,
		numeric_vector & interactionForce, const double rho) {

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
template<size_t dim>
void BGKPseudopotential<dim>::collideAll(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const {

	assert(getForceType() != NO_FORCING);

	if (Stencil_D2Q9 == getStencil()->getStencilType()) {
		collideAllD2Q9(f, densities, velocities, inInitializationProcedure);
	} else {
		// Inefficient collision for other than D2Q9
		BGK::collideAll(f, densities, velocities, locally_owned_dofs,
				inInitializationProcedure);
	}
}

template<size_t dim>
void BGKPseudopotential<dim>::collideAllD2Q9(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		bool inInitializationProcedure) const {

// Efficient collision for D2Q9
	size_t Q = 9;
	size_t D = 2;
	double scaling = getStencil()->getScaling();
	double cs2 = getStencil()->getSpeedOfSoundSquare();
	double prefactor = scaling / cs2;
	double relax_factor = getPrefactor();
	double dt = getDt();

	assert(f.size() == Q);
	assert(velocities.size() == D);
	assert(dim == 2);

	// External force information

	ForceType force_type = getForceType();
	double ext_force_x = getForceX();
	double ext_force_y = getForceY();

#ifdef DEBUG
	size_t n_dofs = f.at(0).size();
	for (size_t i = 0; i < Q; i++) {
		assert (f.at(i).size() == n_dofs);
	}
	assert (densities.size() == n_dofs);
	for (size_t i = 0; i < D; i++) {
		assert (velocities.at(i).size() == n_dofs);
	}
#endif

// allocation
	double feq[9] = { };
	double scalar_product;
	double uSquareTerm;
	double mixedTerm;
	double weighting;
	double rho_i;
	double u_0_i;
	double u_1_i;
	double f_i[9] = { };

	distributed_vector& f0 = f.at(0);
	distributed_vector& f1 = f.at(1);
	distributed_vector& f2 = f.at(2);
	distributed_vector& f3 = f.at(3);
	distributed_vector& f4 = f.at(4);
	distributed_vector& f5 = f.at(5);
	distributed_vector& f6 = f.at(6);
	distributed_vector& f7 = f.at(7);
	distributed_vector& f8 = f.at(8);

	const dealii::DoFHandler<dim>& dof_handler =
			*m_advectionOperator->getDoFHandler();
	const size_t dofs_per_cell = m_advectionOperator->getNumberOfDoFsPerCell();
	std::vector<dealii::types::global_dof_index> localDoFIndices(dofs_per_cell);
	const dealii::UpdateFlags updateFlags = dealii::update_gradients
			| dealii::update_quadrature_points;

	dealii::FEValues<dim> feValues(m_advectionOperator->getMapping(),
			*m_advectionOperator->getFe(),
			*m_advectionOperator->getQuadrature(), updateFlags);

	std::vector<dealii::types::global_dof_index> neighbor_indices(
			dofs_per_cell);
	dealii::FEValues<dim> neighbor_fe_values(m_advectionOperator->getMapping(),
			*m_advectionOperator->getFe(),
			*m_advectionOperator->getQuadrature(), updateFlags);

	densities = 0;
	for (size_t i = 0; i < 9; i++) {
		densities.add(f.at(i));
	}
// loop over all cells
	dealii::Tensor<1, dim> density_gradient;
	typename dealii::DoFHandler<dim>::active_cell_iterator cell =
			dof_handler.begin_active(), endc = dof_handler.end();
	for (; cell != endc; ++cell) {

		if (!cell->is_locally_owned())
			continue;

		// get global degrees of freedom
		cell->get_dof_indices(localDoFIndices);

		// calculate gradients of shape functions
		feValues.reinit(cell);

		// for all dofs in cell
		for (size_t j = 0; j < dofs_per_cell; j++) {

			size_t i = localDoFIndices.at(j);

			// copy fi
			f_i[0] = f0(i);
			f_i[1] = f1(i);
			f_i[2] = f2(i);
			f_i[3] = f3(i);
			f_i[4] = f4(i);
			f_i[5] = f5(i);
			f_i[6] = f6(i);
			f_i[7] = f7(i);
			f_i[8] = f8(i);

			// calculate density
			rho_i = densities(i);//f_i[0] + f_i[1] + f_i[2] + f_i[3] + f_i[4] + f_i[5] + f_i[6]
			//+ f_i[7] + f_i[8];
			//densities(i) = rho_i;
			if (rho_i < 1e-10) {
				throw CollisionException(
						"Densities too small (< 1e-10) for collisions. Decrease time step size.");
			}

			// calculate density gradient
			// This function assumes that the DoFs are uniquely identified with quadrature points
			// it uses the relation grad(rho)(x_i) = grad( sum_alpha f_alpha(x_i))  = grad ( sum_alpha(sum_k(dof_{alpha,k} phi_k(x_i))))
			//                                     = sum_k ((sum_alpha dof_{alpha,k}) grad(phi_k(x_i))) = sum_k (rho_k grad(phi_k(x_i)))

			// PRELIMINATRY IMPLEMENTATION.
			// THE DENSITY GRADIENT CALCULATE MUST BE VERY INEFFICIENT. IT DOES NOT SUPPORT BOUNDARY CONDITIONS AT ALL.
			// IN 3D, ANOTHER LAYER OF CELLS AROUND MUST BE TAKEN CARE OF TO EVALUATE THE GRADIENTS AT VERTICES
			// ONLY GLOBAL REFINEMENT IS SUPPORTED SO FAR
			// WHEN USING A NORMAL FINITE ELEMENT RATHER THAN A DG FINITE ELEMENT, THE SEARCH FOR p ON NEIGHBOR ELEMENT COULD
			// BE SPEEDED UP

			// average face gradients over cell and all neighbors that also contain the support point
			density_gradient = 0;
/*			double weight_sum = 0;
			dealii::Point<dim> p = feValues.quadrature_point(
					celldofToQIndex.at(j));
			// a list of all cells that contain p
			vector<typename dealii::DoFHandler<dim>::cell_iterator> cells_with_point;
			// this cell
			cells_with_point.push_back(cell);
			// iterate over neighbors
			for (size_t k = 0; k < dealii::GeometryInfo<dim>::faces_per_cell;
					k++) {
				if (cell->face(k)->at_boundary()) {
					// TODO incorporate periodic boundaries
					continue;
				}
				if (cell->neighbor(k)->point_inside(p)) {
					cells_with_point.push_back(cell->neighbor(k));
				}
			}
			// iterate over neighbors of the neighbors that are in the list
			size_t len = cells_with_point.size();
			for (size_t c = 1; c < len; c++) {
				for (size_t k = 0;
						k < dealii::GeometryInfo<dim>::faces_per_cell; k++) {
					if (cells_with_point.at(c)->face(k)->at_boundary()) {
						// TODO incorporate periodic boundaries
						continue;
					}
					bool already_visited = (std::find(cells_with_point.begin(),
							cells_with_point.end(),
							cells_with_point.at(c)->neighbor(k))
							!= cells_with_point.end());
					if ((cells_with_point.at(c)->neighbor(k)->point_inside(p))
							and (not already_visited)) {
						cells_with_point.push_back(
								cells_with_point.at(c)->neighbor(k));
					}
				}

			}
			// calculate the density gradient as the average over all cells that contain p
			//size_t n_finds = 0;
			len = cells_with_point.size();
			for (size_t c = 0; c < len; c++) {
				neighbor_fe_values.reinit(cells_with_point.at(c));
				cells_with_point.at(c)->get_dof_indices(neighbor_indices);

				 // find point in neighbor quadrature
				 size_t m;
				 bool found = false;
				 for (m = 0; m < feValues.dofs_per_cell; m++) {
				 if (p.distance(
				 neighbor_fe_values.quadrature_point(
				 celldofToQIndex.at(m))) < 1e-12) {
				 n_finds++;
				 found = true;
				 break;
				 }
				 }
				 assert(found);// otherwise: the point is not a quadrature point at neighbor cell.
				 // requires more complicated evaluation.
				 *

				// calculate average gradient
				for (size_t n = 0; n < dofs_per_cell; n++) {
					double rho_n = densities(neighbor_indices.at(n));
					for (size_t m = 0; m < dofs_per_cell; m++) {
						double weight = exp(-p.distance(
								neighbor_fe_values.quadrature_point(
										celldofToQIndex.at(m)))
								/ cell->diameter());
						density_gradient += (neighbor_fe_values.shape_grad(n,
								celldofToQIndex.at(m)) * rho_n * weight);
						weight_sum += weight;
					}
				}
			}
			cout << endl;
			// average
			density_gradient *= (1. / (weight_sum));
			double limit = 5;
			 if (density_gradient.norm() > limit){
			 density_gradient*=(limit/density_gradient.norm());
			 }
			 */
			/*cout  << rho_i << " | "
					<< density_gradient[0] << " " << density_gradient[1]
					<< endl;
*/
			// =============================
			// Calculate Interaction Force
			// =============================

			// Determine constant force
//			double forceX = 16e-5;			//7.92*(1e-6) ;
///			double forceY = (double) 0;

			double force_x = 0;
			double force_y = 0;
			if (SHAN_CHEN == m_pseudopotentialParameters.pseudopotentialType) {
				// Psi = 1-exp(-rho)
				double psi = (double) 1 - exp(-rho_i);
				double G = m_pseudopotentialParameters.G;
				force_x = -G * psi * exp(-rho_i) * density_gradient[0]
						+ ext_force_x;
				force_y = -G * psi * exp(-rho_i) * density_gradient[1]
						+ ext_force_y;
			} else if (SUKOP
					== m_pseudopotentialParameters.pseudopotentialType) {
				// Psi = psi0 * exp(-rho/rho0)
				double G = m_pseudopotentialParameters.G;
				double psi0 = (double) 4;
				double rho0 = (double) 200;
				double psi = (double) psi0 * exp(-rho_i / rho0);
				force_x = -G * psi * (-psi0 / rho0) * exp(-rho_i / rho0)
						* density_gradient[0] + ext_force_x;
				force_y = -G * psi * (-psi0 / rho0) * exp(-rho_i / rho0)
						* density_gradient[1] + ext_force_y;
			} else if (CARNAHAN_STARLING
					== m_pseudopotentialParameters.pseudopotentialType) {
				// Psi Carnahan Starling
				double G = m_pseudopotentialParameters.G;
				double T = m_pseudopotentialParameters.T * cs2;
				double pCS = rho_i * T
						* (1. + rho_i + pow(rho_i, 2.) - pow(rho_i, 3.))
						/ pow(1. - rho_i, 3.) - pow(rho_i, 2.);
				double psi = sqrt(2. * (pCS - cs2 * rho_i) / G);
				double grad_pCS_X = density_gradient[0]
						* ((T
								* (1 + 4 * rho_i + 4 * pow(rho_i, 2.)
										- 4 * pow(rho_i, 3.) + pow(rho_i, 4.))
								/ pow(1. - rho_i, 4.)) - (2. * rho_i));
				double grad_pCS_Y = density_gradient[1]
						* ((T
								* (1 + 4 * rho_i + 4 * pow(rho_i, 2.)
										- 4 * pow(rho_i, 3.) + pow(rho_i, 4.))
								/ pow(1. - rho_i, 4.)) - (2. * rho_i));
				double gradPsiX = grad_pCS_X
						- cs2 * density_gradient[0]
								/ (2. * G * (pCS - cs2 * rho_i));
				double gradPsiY = grad_pCS_Y
						- cs2 * density_gradient[1]
								/ (2. * G * (pCS - cs2 * rho_i));
				force_x = -G * psi * gradPsiX + ext_force_x;
				force_y = -G * psi * gradPsiY + ext_force_y;
			}
			// With Tuneable surface tension (Kupershtokh)
//			double A = 0.0 ;
//			double psi = 1-exp(-rho_i) ;
//			double force_x = -m_G*(2*A*exp(-rho_i)*density_gradient[0] + (1-A)*psi*exp(-rho_i)*density_gradient[0]) + ext_force_x;
//			double force_y = -m_G*(2*A*exp(-rho_i)*density_gradient[1] + (1-A)*psi*exp(-rho_i)*density_gradient[1]) + ext_force_y;
//			double force_x = -m_G * 0.5 * density_gradient[0] * exp(-rho_i) * (1-exp(-rho_i)) * (A+1) + ext_force_x;
//			double force_y = -m_G * 0.5 * density_gradient[1] * exp(-rho_i) * (1-exp(-rho_i)) * (A+1) + ext_force_y;

			if (not inInitializationProcedure) {
				// calculate macroscopic velocity (velocities.at()(i)) and equilibrium velocity (u_0_i, u_1_i)
				// for all velocity components
				u_0_i = scaling / rho_i
						* (f_i[1] + f_i[5] + f_i[8] - f_i[3] - f_i[6] - f_i[7]);
				u_1_i = scaling / rho_i
						* (f_i[2] + f_i[5] + f_i[6] - f_i[4] - f_i[7] - f_i[8]);
				if (force_type == NO_FORCING) {
					velocities.at(0)(i) = u_0_i;
					velocities.at(1)(i) = u_1_i;
				} else if (force_type == SHIFTING_VELOCITY) {
					velocities.at(0)(i) = u_0_i + 0.5 * dt * force_x / rho_i;
					velocities.at(1)(i) = u_1_i + 0.5 * dt * force_y / rho_i;
					u_0_i -= (double) 1 / relax_factor * dt * force_x / rho_i;
					u_1_i -= (double) 1 / relax_factor * dt * force_y / rho_i;
				} else if (force_type == EXACT_DIFFERENCE) {
					velocities.at(0)(i) = u_0_i + 0.5 * dt * force_x / rho_i;
					velocities.at(1)(i) = u_1_i + 0.5 * dt * force_y / rho_i;
				} else { // GUO
					u_0_i += 0.5 * dt * force_x / rho_i;
					u_1_i += 0.5 * dt * force_y / rho_i;
					velocities.at(0)(i) = u_0_i;
					velocities.at(1)(i) = u_1_i;
				}
			}

			// Calculate equilibrium distribution
			scalar_product = u_0_i * u_0_i + u_1_i * u_1_i;
			uSquareTerm = -scalar_product / (2 * cs2);
			// direction 0
			weighting = 4. / 9. * rho_i;
			feq[0] = weighting * (1 + uSquareTerm);
			// directions 1-4
			weighting = 1. / 9. * rho_i;
			mixedTerm = prefactor * (u_0_i);
			feq[1] = weighting
					* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
			feq[3] = weighting
					* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
			mixedTerm = prefactor * (u_1_i);
			feq[2] = weighting
					* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
			feq[4] = weighting
					* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
			// directions 5-8
			weighting = 1. / 36. * rho_i;
			mixedTerm = prefactor * (u_0_i + u_1_i);
			feq[5] = weighting
					* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
			feq[7] = weighting
					* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
			mixedTerm = prefactor * (-u_0_i + u_1_i);
			feq[6] = weighting
					* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
			feq[8] = weighting
					* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);

			// BGK collision
			f_i[0] += relax_factor * (f_i[0] - feq[0]);
			f_i[1] += relax_factor * (f_i[1] - feq[1]);
			f_i[2] += relax_factor * (f_i[2] - feq[2]);
			f_i[3] += relax_factor * (f_i[3] - feq[3]);
			f_i[4] += relax_factor * (f_i[4] - feq[4]);
			f_i[5] += relax_factor * (f_i[5] - feq[5]);
			f_i[6] += relax_factor * (f_i[6] - feq[6]);
			f_i[7] += relax_factor * (f_i[7] - feq[7]);
			f_i[8] += relax_factor * (f_i[8] - feq[8]);

			// Add Source term
			// Exact difference method (Kupershtokh)
			if (force_type == EXACT_DIFFERENCE) {
				ExternalForceFunctions::applyExactDifferenceForcingD2Q9(f_i,
						force_x, force_y, u_0_i, u_1_i, rho_i, getDt(),
						prefactor);
			}
			// TODO add source term guo

			// copy to global variables
			f0(i) = f_i[0];
			f1(i) = f_i[1];
			f2(i) = f_i[2];
			f3(i) = f_i[3];
			f4(i) = f_i[4];
			f5(i) = f_i[5];
			f6(i) = f_i[6];
			f7(i) = f_i[7];
			f8(i) = f_i[8];

		} /* for all dofs */
	} /* for all cells */
} /* collideAllD2Q9 */

// Explicit Instantiation
template class BGKPseudopotential<2> ;
template class BGKPseudopotential<3> ;
} /* namespace natrium */
