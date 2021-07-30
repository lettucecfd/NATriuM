/*
 * SemiLagrangianBoundaryHandler.h
 *
 *  Created on: 15.11.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_BOUNDARIES_SEMILAGRANGIANBOUNDARYHANDLER_H_
#define LIBRARY_NATRIUM_BOUNDARIES_SEMILAGRANGIANBOUNDARYHANDLER_H_

#include <map>
#include <utility>

#include "deal.II/dofs/dof_handler.h"
#include "deal.II/fe/fe_update_flags.h"
#include "deal.II/fe/fe_values.h"
#include "deal.II/base/point.h"

#include "../advection/SemiLagrangianTools.h"
#include "../problemdescription/BoundaryCollection.h"
#include "../utilities/BasicNames.h"
#include "BoundaryHit.h"

namespace natrium {

// forward declaration
template<size_t dim>
class AdvectionOperator;

template<size_t dim>
class SemiLagrangian;

/**
 * @short Handler for semi-Lagrangian boundary conditions. This class manages the update of distribution functions
 *		 at the boundaries in the semi-Lagrangian streaming step.
 */
template<size_t dim>
class SemiLagrangianBoundaryHandler {
private:

	/// list of boundary hit points
	HitList<dim> m_hitList;

	/// dt
	double m_timeStep;

	/// the LBM stencil (e.g. D2Q9)
	const Stencil& m_stencil;

	/// the advection operator (usually a semi-Lagrangian one)
	const AdvectionOperator<dim>& m_advection;

	/// definition of the flow
	const ProblemDescription<dim>& m_problem;

	/// definition of the boundary conditions
	const BoundaryCollection<dim>& m_boundaries;

	/// deal.II's degree of freedom handler
	const dealii::DoFHandler<dim>& m_dof;

public:
	/**
	 * @short
	 */
	SemiLagrangianBoundaryHandler(SemiLagrangian<dim>& advection) :
			m_timeStep(advection.getDeltaT()), m_stencil(
					*advection.getStencil()), m_advection(advection), m_problem(
					advection.getProblem()), m_boundaries(
					*advection.getBoundaries()), m_dof(
					*advection.getDoFHandler()) {
		assert(m_stencil.getD() == dim);
		assert(m_timeStep >= 0);

	}
	virtual ~SemiLagrangianBoundaryHandler() {

	}
	void addHit(const LagrangianPathTracker<dim>& tracker, size_t boundary_id);

	void apply(DistributionFunctions& f_new, const DistributionFunctions& f_old, double t);

	void applyToG(DistributionFunctions& f, DistributionFunctions& g, double t, const double gamma);

	/*void print_out(){
	 pout << "Semi Lagrangian Boundary Handler Object:" << endl;
	 m_hitList.print_out();
	 pout << "----------------------------------------" << endl;
	 pout << "#Hits at support points: " << endl;
	 }*/

	double getTimeStep() const {
		return m_timeStep;
	}

	void setTimeStep(double timeStep) {
		m_timeStep = timeStep;
	}

	size_t n_hits() const {
		return m_hitList.n_hits();
	}

	size_t n_cells() const {
		return m_hitList.n_cells();
	}
}
;

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_BOUNDARIES_SEMILAGRANGIANBOUNDARYHANDLER_H_ */
