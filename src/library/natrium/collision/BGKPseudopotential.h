/*
 * BGKPseudopotential.h
 *
 *  Created on: Nov 11, 2014
 *      Author: kk
 */

#ifndef BGKPSEUDOPOTENTIAL_H_
#define BGKPSEUDOPOTENTIAL_H_

#include "BGK.h"

#include "deal.II/grid/tria.h"
#include "deal.II/numerics/vector_tools.h"

//#include "../advection/SEDGMinLee.h"

#include "../utilities/BasicNames.h"

namespace natrium {

// forward declaration
template<size_t dim>
class SEDGMinLee;
class Stencil;

class BGKPseudopotential: public BGK {
private:
	shared_ptr<SEDGMinLee<2> > m_advectionOperator;

public:

	/// constructor
	BGKPseudopotential(double relaxationParameter, double dt, const shared_ptr<Stencil> stencil);

	/// destructor
	virtual ~BGKPseudopotential();

	virtual double getEquilibriumDistribution(size_t i, const numeric_vector& u,
			const double rho = 1) const;

	/**
	 * @short function for collision
	 * @short f the global vectors of discrete particle distribution functions
	 * @short densities the global vector of densities
	 * @short velocities the global vectors of velocity components [ [u_1x, u_2x, ...], [u_1y, u_2y, ...] ]
	 * @short inInitializationProcedure indicates if the collision is performed in the context of an iterative initilizatation procedure. In this case, only the macroscopic densities are recalculated, while the velocities remain unchanged. default: false
	 */
	virtual void collideAll(DistributionFunctions& f,
			distributed_vector& densities,
			vector<distributed_vector>& velocities,
			bool inInitializationProcedure = false,
			const dealii::IndexSet& locally_owned_dofs = dealii::IndexSet()) const;

	void getInteractionForce(const vector<double>& distributions,
			numeric_vector & interactionForce, const double rho = 1);

	const shared_ptr<SEDGMinLee<2> > getAdvectionOperator() const {
		return m_advectionOperator;
	}

	void setAdvectionOperator(
			const shared_ptr<SEDGMinLee<2> >& advectionOperator) {
		m_advectionOperator = advectionOperator;
	}
};

} /* namespace natrium */

#endif /* BGKPSEUDOPOTENTIAL_H_ */
