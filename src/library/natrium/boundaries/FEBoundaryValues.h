/*
 * FEBoundaryValues.h
 *
 *  Created on: 29.11.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_BOUNDARIES_FEBOUNDARYVALUES_H_
#define LIBRARY_NATRIUM_BOUNDARIES_FEBOUNDARYVALUES_H_

#include "deal.II/fe/fe_values.h"

#include "BoundaryFlags.h"
#include "BoundaryTools.h"
#include "../utilities/BasicNames.h"
#include "../advection/SemiLagrangianTools.h"

namespace natrium {

/**
 * @short A class that calculates macroscopic values at the boundaries,
 * 			implemented similarly as deal.II's FEValues class.
 * 			In constrast to the latter, a new instance of FEBoundaryValues
 * 			is to be created for each cell by the semi-Lagrangian boundary
 * 			handler. This is due to the fact that each cell may have
 * 			a new set of boundary hit points. The instance is initialized for
 * 			the entire cell. Boundary values at single quadrature points
 * 			are obtained by calling reinit(q_point) for the point and then
 * 			use the getter-Functions to retrieve the desired values.
 */
template<size_t dim>
class FEBoundaryValues {
private:

	// underlying fe values object
	dealii::FEValues<dim> m_feValues;
	// boundary flags of all values that are to be calculated from the boundary distribution functions
	BoundaryFlags m_flags;
	// all flags that are up to date after calling reinit(), mainly for debugging purposes
	BoundaryFlags m_upToDate;

	const GlobalBoundaryData& m_data;

	// constants
	size_t m_dofsPerCell;
	std::vector<dealii::types::global_dof_index> m_localDofs;

	// macroscopic variables and dof variables
	double m_rho;
	dealii::Tensor<1, dim> m_rho_u;
	dealii::Tensor<1, dim> m_u;
	vector<double> m_rho_dof;
	vector<double> m_distributions;
	vector<dealii::Tensor<1, dim> > m_rho_u_dof;
	vector<dealii::Tensor<1, dim> > m_u_dof;

	// wall values and temporal derivatives
	//double rho_w;
	double m_drhodt;
	//dealii::Tensor<1, dim> u_w;
	dealii::Tensor<1, dim> m_dudt;
	//dealii::Tensor<2, dim> dudx_w;
	//dealii::Tensor<2, dim> d2udxdt_w;
public:
	/**
	 * @short Constructor
	 * @param cell the current cell. Note that every cell requires a new instance of FEBoundaryValues for the update.
	 * @param hit_points the local boundary hit points
	 * @param fe the (Lagrangian) finite elements that are used for discretization
	 * @param mapping mapping from real to unit cell
	 * @param flags the boundary update flags
	 * @param global_data GlobalData instance that gives quick access to relevant data, such as the time step, the stencil and the distributions
	 */
	FEBoundaryValues(
			const typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
			const vector<dealii::Point<dim> >& hit_points,
			const dealii::FiniteElement<dim>& fe,
			const dealii::Mapping<dim>& mapping, BoundaryFlags& flags,
			const GlobalBoundaryData& global_data) :
			m_feValues(mapping, fe,
					makeQuadratureAtPoints<dim>(cell, hit_points, mapping),
					getFEValuesFlags(flags)), m_flags(flags), m_data(
					global_data) {
		m_dofsPerCell = fe.n_dofs_per_cell();
		resize(m_dofsPerCell);
		setToZero();
		cell->get_dof_indices(m_localDofs);
	}


	virtual ~FEBoundaryValues(){

	}

	/**
	 * @short resize all vectors
	 */
	void resize(size_t n) {
		m_rho_dof.resize(n);
		m_rho_u_dof.resize(n);
		m_u_dof.resize(n);
		m_localDofs.resize(n);
	}

	void setToZero() {
		// macroscopic variables and dof variables
		m_rho = 0;
		m_rho_u = 0;
		m_u = 0;

		// wall values and temporal derivatives
		//rho_w = 0;
		m_drhodt = 0;
		//u_w = 0;
		m_dudt = 0;
		//dudx_w = 0;
		//d2udxdt_w = 0;
	}

	void clear() {
		m_dofsPerCell = 0;

		// macroscopic variables and dof variables
		m_rho = 0;
		m_rho_u = 0;
		m_u = 0;
		m_rho_dof.clear();
		m_rho_u_dof.clear();
		m_u_dof.clear();
		m_localDofs.clear();

		// wall values and temporal derivatives
		//rho_w = 0;
		m_drhodt = 0;
		//u_w = 0;
		m_dudt = 0;
		//dudx_w = 0;
		//d2udxdt_w = 0;
	}

	void reinit(size_t q_point) {
		m_upToDate = only_distributions;

		// Has to be in the right order:
		// calculate rho
		if (m_flags != only_distributions) {
			calculateRho(q_point);
		}
		// calculate u
		if ((boundary_drho_dt & m_flags) or (boundary_u & m_flags)
				or (boundary_du_dt & m_flags)) {
			calculateU(q_point);
		}
		// calculate drho/dt
		if (boundary_drho_dt & m_flags) {
			calculateDRhoDt(q_point);
		}
		// calculate du/dt
		if (boundary_du_dt & m_flags) {
			calculateDuDt(q_point);
		}
		assert(m_flags == m_upToDate);
	}

	size_t getDofsPerCell() const {
		return m_dofsPerCell;
	}

	double getDrhodt() const {
		assert(boundary_drho_dt & m_upToDate);
		return m_drhodt;
	}

	const dealii::Tensor<1, dim>& getDudt() const {
		return m_dudt;
	}

	BoundaryFlags getFlags() const {
		return m_flags;
	}

	const std::vector<dealii::types::global_dof_index>& getLocalDofs() const {
		return m_localDofs;
	}

	double getRho() const {
		return m_rho;
	}

	const vector<double>& getRhoDof() const {
		return m_rho_dof;
	}

	const dealii::Tensor<1, dim>& getRhoU() const {
		return m_rho_u;
	}

	const vector<dealii::Tensor<1, dim> >& getRhoUDof() const {
		return m_rho_u_dof;
	}

	const dealii::Tensor<1, dim>& getU() const {
		return m_u;
	}

	const vector<dealii::Tensor<1, dim> >& getUDof() const {
		return m_u_dof;
	}

	BoundaryFlags getUpToDate() const {
		return m_upToDate;
	}

	double getDistribution(size_t i, size_t q_point) {
		assert(m_feValues.get_fe().n_dofs_per_cell() == m_dofsPerCell);
		assert(m_dofsPerCell == m_localDofs.size());
		assert(i < m_data.m_Q);
		// calculate distribution i at q_point and (t-delta_t)
		double result = 0;
		for (size_t dof = 0; dof < m_dofsPerCell; dof++) {
			result += m_data.m_fold.at(i)(m_localDofs.at(dof))
					* m_feValues.shape_value(dof, q_point);
		}
		return result;
	}

	const dealii::Point<dim>& getPoint(size_t q_point){
		return m_feValues.get_quadrature().point(q_point);
	}

	const GlobalBoundaryData& getData() const {
		return m_data;
	}

private:

	dealii::UpdateFlags getFEValuesFlags(BoundaryFlags& flags) {
		dealii::UpdateFlags result = dealii::update_values;
		if (boundary_drho_dt & flags) {
			result |= dealii::update_gradients;
		}
		if (boundary_du_dt & flags) {
			result |= dealii::update_gradients;
			result |= dealii::update_hessians;
		}
		return result;
	}
	/**
	 * @short calculate boundary density at last time step (t - delta_t)
	 */
	void calculateRho(size_t q_point) {

		// check that object has been initialized properly via reinit
		assert(not (boundary_rho & m_upToDate));
		assert(m_data.m_viscosity > 0);
		assert(m_data.m_dt > 0); // time step

		// check function parameters
		assert(m_feValues.get_fe().n_dofs_per_cell() == m_dofsPerCell);
		assert(m_dofsPerCell == m_localDofs.size());

		// initalize macroscopic variables and dof variables
		m_rho = 0;
		for (size_t i = 0; i < m_dofsPerCell; i++) {
			m_rho_dof.at(i) = 0;
		}

		// calculate macroscopic variables and dof variables
		for (size_t dof = 0; dof < m_dofsPerCell; dof++) {
			for (size_t qi = 0; qi < m_data.m_Q; qi++) {
				m_rho_dof.at(dof) += m_data.m_fold.at(qi)(m_localDofs.at(dof));
			}
			m_rho += m_rho_dof.at(dof) * m_feValues.shape_value(dof, q_point);
		}
		m_upToDate |= boundary_rho;
	}

	/**
	 * @short calculate boundary velocity at last time step (t - delta_t)
	 */
	void calculateU(size_t q_point) {
		// check that density is already updated
		assert(boundary_rho & m_upToDate);
		// check that velocity is not updated
		assert(not (boundary_u & m_upToDate));
		// check input
		assert(m_data.m_viscosity > 0);
		assert(m_data.m_dt > 0); // time step

		// check function parameters
		assert(m_feValues.get_fe().n_dofs_per_cell() == m_dofsPerCell);
		assert(m_dofsPerCell == m_localDofs.size());

		// initalize macroscopic variables and dof variables
		m_rho_u = 0;
		m_u = 0;
		for (size_t i = 0; i < m_dofsPerCell; i++) {
			m_rho_u_dof.at(i) = 0;
			m_u_dof.at(i) = 0;
		}

		// calculate macroscopic variables and dof variables
		for (size_t dof = 0; dof < m_dofsPerCell; dof++) {
			for (size_t qi = 0; qi < m_data.m_Q; qi++) {
				for (size_t i = 0; i < dim; i++) {
					m_rho_u_dof.at(dof)[i] += m_data.m_fold.at(qi)(
							m_localDofs.at(dof))
							* m_data.m_stencil.getDirection(qi)[i];
				}
			}
			m_u_dof.at(dof) = m_rho_u_dof.at(dof);
			m_u_dof.at(dof) /= m_rho_dof.at(dof);
			m_rho_u += m_rho_u_dof.at(dof)
					* m_feValues.shape_value(dof, q_point);
			m_u += m_u_dof.at(dof) * m_feValues.shape_value(dof, q_point);
		}
		assert((m_u * m_rho - m_rho_u).norm() < 1e-10);
		m_upToDate |= boundary_u;
	}

	/**
	 * @short calculate temporal derivative of rho at last time step (t-delta_t) by compressible conti equation
	 */
	void calculateDRhoDt(size_t q_point) {

		assert(not (boundary_drho_dt & m_upToDate));
		assert(boundary_u & m_upToDate);

		m_drhodt = 0;

		// calculate divergence
		for (size_t dof = 0; dof < m_dofsPerCell; dof++) {
			for (size_t i = 0; i < dim; i++) {
				// derivative drhoui_dxi
				m_drhodt += m_feValues.shape_grad(dof, q_point)[i]
						* m_rho_u_dof.at(dof)[i];
			}
		}
		m_drhodt *= -1;
		m_upToDate |= boundary_drho_dt;
	}

	/**
	 * @short calculate temporal derivative of u at last time step (t-delta_t) by compressible momentum equation
	 */
	void calculateDuDt(size_t q_point) {

		assert(not (boundary_du_dt & m_upToDate));
		assert(boundary_u & m_upToDate);
		assert(boundary_rho & m_upToDate);

		m_dudt = 0;

		// calculate u_w: first viscous term, then derivative part, then rest
		for (size_t i = 0; i < dim; i++) {
			// viscous term
			for (size_t dof = 0; dof < m_dofsPerCell; dof++) {
				for (size_t k = 0; k < dim; k++) {
					m_dudt[i] += m_feValues.shape_hessian(dof, q_point)[k][k]
							* m_u_dof.at(dof)[i];
					m_dudt[i] += m_feValues.shape_hessian(dof, q_point)[i][k]
							* m_u_dof.at(dof)[k];
				}
			}
			m_dudt[i] *= (-m_data.m_viscosity);
			for (size_t dof = 0; dof < m_dofsPerCell; dof++) {
				// nonlinear term
				for (size_t k = 0; k < dim; k++) {
					m_dudt[i] += m_u[k] * m_feValues.shape_grad(dof, q_point)[k]
							* m_u_dof.at(dof)[i];
				}
				// pressure term
				m_dudt[i] += m_data.m_cs2 / m_rho
						* m_feValues.shape_grad(dof, q_point)[i]
						* m_rho_dof.at(dof);
			}
			m_dudt[i] *= -1;
		}

		m_upToDate |= boundary_du_dt;
	}

	// TODO (possibly) calculate du_dx d2u_dxdt
	/**
	 * @short calculate temporal derivative of rho at last time step (t-delta_t) by compressible conti equation
	 */
	/*
	 template<size_t dim>
	 void calculateWallValues(const GlobalBoundaryData& g, size_t q_point) {

	 assert(not (boundary_drho_dt & m_upToDate));
	 assert(boundary_u & m_upToDate);

	 b.rho_w = 0;
	 b.drhodt_w = 0;
	 b.u_w = 0;
	 b.dudt_w = 0;
	 b.dudx_w = 0;
	 b.d2udxdt_w = 0;

	 // calculate rho_w: first divergence, then rest
	 for (size_t dof = 0; dof < b.dofs_per_cell; dof++) {
	 for (size_t i = 0; i < dim; i++) {
	 // derivative drhoui_dxi
	 b.drhodt_w += fe_values.shape_grad(dof, q_point)[i]
	 * b.rho_u_dof.at(dof)[i];
	 }
	 }
	 b.drhodt_w *= -1;

	 // calculate u_w: first viscous term, then derivative part, then rest
	 for (size_t i = 0; i < dim; i++) {
	 // viscous term
	 for (size_t dof = 0; dof < b.dofs_per_cell; dof++) {
	 for (size_t k = 0; k < dim; k++) {
	 b.dudt_w[i] += fe_values.shape_hessian(dof, q_point)[k][k]
	 * b.u_dof.at(dof)[i];
	 b.dudt_w[i] += fe_values.shape_hessian(dof, q_point)[i][k]
	 * b.u_dof.at(dof)[k];
	 }
	 }
	 b.dudt_w[i] *= (-g.m_viscosity);
	 for (size_t dof = 0; dof < b.dofs_per_cell; dof++) {
	 // nonlinear term
	 for (size_t k = 0; k < dim; k++) {
	 b.dudt_w[i] += b.u[k] * fe_values.shape_grad(dof, q_point)[k]
	 * b.u_dof.at(dof)[i];
	 }
	 // pressure term
	 b.dudt_w[i] += g.m_cs2 / b.rho * fe_values.shape_grad(dof, q_point)[i]
	 * b.rho_dof.at(dof);
	 }
	 b.u_w[i] *= -1;
	 }

	 // calculate dudx_w
	 for (size_t i = 0; i < dim; i++) {
	 for (size_t j = 0; j < dim; j++) {
	 // viscous term
	 for (size_t dof = 0; dof < b.dofs_per_cell; dof++) {
	 for (size_t k = 0; k < dim; k++) {
	 b.d2udxdt_w[i][j] += fe_values.shape_3rd_derivative(dof,
	 q_point)[k][k][j] * b.u_dof.at(dof)[i];
	 b.d2udxdt_w[i][j] += fe_values.shape_3rd_derivative(dof,
	 q_point)[i][k][j] * b.u_dof.at(dof)[k];
	 }
	 }
	 b.d2udxdt_w[i][j] *= (-g.m_viscosity);
	 for (size_t dof = 0; dof < b.dofs_per_cell; dof++) {
	 for (size_t k = 0; k < dim; k++) {
	 // nonlinear term 1
	 b.d2udxdt_w[i][j] += fe_values.shape_grad(dof, q_point)[j]
	 * b.u_dof.at(dof)[k]
	 * fe_values.shape_grad(dof, q_point)[k]
	 * b.u_dof.at(dof)[i];
	 // nonlinear term 2
	 b.d2udxdt_w[i][j] += b.u[j]
	 * fe_values.shape_hessian(dof, q_point)[j][k]
	 * b.u_dof.at(dof)[i];

	 }
	 // pressure term 1
	 b.d2udxdt_w[i][j] += g.m_cs2 / (b.rho * b.rho)
	 * fe_values.shape_hessian(dof, q_point)[i][j]
	 * b.rho_dof.at(dof);
	 // pressure term 2
	 b.d2udxdt_w[i][j] += g.m_cs2 / (b.rho * b.rho)
	 * fe_values.shape_grad(dof, q_point)[i]
	 * b.rho_dof.at(dof)
	 * fe_values.shape_grad(dof, q_point)[j]
	 * b.rho_dof.at(dof);
	 }
	 b.d2udxdt_w[i] *= -1;
	 // + dudx_w, which is not stored, but calculated by:
	 for (size_t dof = 0; dof < b.dofs_per_cell; dof++) {
	 b.dudx_w[i][j] += fe_values.shape_grad(dof, q_point)[j]
	 * b.u_dof.at(dof)[i];
	 }
	 }
	 }
	 }
	 */

};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_BOUNDARIES_FEBOUNDARYVALUES_H_ */
