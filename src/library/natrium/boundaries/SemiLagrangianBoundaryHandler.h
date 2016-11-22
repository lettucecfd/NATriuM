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


/**
 * @short Handler for semi-Lagrangian boundary conditions.
 */
template<size_t dim>
class SemiLagrangianBoundaryHandler {
private:
	HitList<dim> m_hitList;
	double m_timeStep;
	const Stencil& m_stencil;
	const BoundaryCollection<dim>& m_boundaries;
public:
	SemiLagrangianBoundaryHandler(double dt, const Stencil& stencil,
			const BoundaryCollection<dim>& boundaries) :
			m_timeStep(dt), m_stencil(stencil), m_boundaries(boundaries) {
		assert(stencil.getD() == dim);
		assert(dt >= 0);

	}
	virtual ~SemiLagrangianBoundaryHandler() {

	}
	void addHit(const LagrangianPathTracker<dim>& tracker, size_t boundary_id,
			const AdvectionOperator<dim>& sl);

	void apply(DistributionFunctions& f_new, const DistributionFunctions& f_old,
			const dealii::DoFHandler<dim>& dof,
			boost::shared_ptr<dealii::FEValues<dim> > fe_values, double t);

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
