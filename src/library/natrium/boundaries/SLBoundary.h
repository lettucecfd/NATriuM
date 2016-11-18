/*
 * SLBoundary.h
 *
 *  Created on: 19.09.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_BOUNDARIES_SLBOUNDARY_H_
#define LIBRARY_NATRIUM_BOUNDARIES_SLBOUNDARY_H_

#include "Boundary.h"

#include "deal.II/base/function.h"

#include "../solver/DistributionFunctions.h"
#include "../stencils/Stencil.h"
#include "../advection/SemiLagrangianTools.h"
#include "../utilities/BasicNames.h"
#include "../advection/AdvectionOperator.h"

namespace natrium {

template<size_t dim>
class SLBoundary: public Boundary<dim> {
private:
	size_t m_boundaryIndicator;

	boost::shared_ptr<dealii::Function<dim> > m_boundaryDensity;

	boost::shared_ptr<dealii::Function<dim> > m_boundaryVelocity;

public:

	SLBoundary(size_t boundaryIndicator) :
			m_boundaryIndicator(boundaryIndicator) {

	}
	// TODO boundary indicator
	SLBoundary(size_t boundaryIndicator,
			boost::shared_ptr<dealii::Function<dim> > boundaryDensity,
			boost::shared_ptr<dealii::Function<dim> > boundaryVelocity) :
			m_boundaryIndicator(boundaryIndicator), m_boundaryDensity(
					boundaryDensity), m_boundaryVelocity(boundaryVelocity) {
	}
	;

	virtual ~SLBoundary() {
	}
	;

	virtual bool isPeriodic() const {
		return false;
	}
	virtual bool isLinearFluxBoundary() const {
		return false;
	}
	virtual bool isSLBoundary() const {
		return true;
	}

	size_t getBoundaryIndicator() const {
		return m_boundaryIndicator;
	}
	const boost::shared_ptr<dealii::Function<dim> >& getBoundaryDensity() const {
		return m_boundaryDensity;
	}

	const boost::shared_ptr<dealii::Function<dim> >& getBoundaryVelocity() const {
		return m_boundaryVelocity;
	}

	virtual void calculateBoundaryValues(const DistributionFunctions& f_old,
			DistributionFunctions& f_new,
			const std::vector< dealii::types::global_dof_index > & local_dofs,
			const dealii::FEValues<dim>&,
			size_t q_point, const LagrangianPathDestination& destination,
			double dt) const = 0;

	virtual dealii::UpdateFlags getUpdateFlags() const = 0;
};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_BOUNDARIES_SLBOUNDARY_H_ */

