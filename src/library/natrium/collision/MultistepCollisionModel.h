/*
 * MultistepCollisionModel.h
 *
 *  Created on: 20.06.2016
 *      Author: dominik
 */

#ifndef MULTISTEPCOLLISIONMODEL_H_
#define MULTISTEPCOLLISIONMODEL_H_

#include <natrium/collision/CollisionModel.h>

namespace natrium {

class MultistepCollisionModel: public CollisionModel {

private:
	DistributionFunctions m_formerF;
	DistributionFunctions m_formerFEq;

public:
	MultistepCollisionModel(double relaxationParameter, double dt, const boost::shared_ptr<Stencil> stencil);


	virtual ~MultistepCollisionModel();

	virtual void collideAll(DistributionFunctions& f,
				distributed_vector& densities,
				vector<distributed_vector>& velocities,
				const dealii::IndexSet& locally_owned_dofs,
				bool inInitializationProcedure = false) const;

	double getEquilibriumDistribution(size_t i, const numeric_vector& u,
				const double rho) const;

	DistributionFunctions getFormerF()
	{
		return m_formerF;
	}

	DistributionFunctions getFormerFEq()
	{
		return m_formerFEq;
	}



};

} /* namespace natrium */

#endif /* MULTISTEPCOLLISIONMODEL_H_ */
