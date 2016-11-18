/*
 * MultistepCollisionData.h
 *
 *  Created on: 20.06.2016
 *      Author: dominik
 */

#ifndef MULTISTEPCOLLISIONDATA_H_
#define MULTISTEPCOLLISIONDATA_H_

#include "../solver/DistributionFunctions.h"

namespace natrium {

class MultistepCollisionData { //: public CollisionModel {

private:
	// double m_dt;


protected:
	mutable DistributionFunctions m_formerF;
	mutable DistributionFunctions m_formerFEq;
	mutable bool m_firstCollision = 1;

public:
	MultistepCollisionData();//double relaxationParameter, double dt, const boost::shared_ptr<Stencil> stencil);


	virtual ~MultistepCollisionData();

/*	virtual void collideAll(DistributionFunctions& f,
				distributed_vector& densities,
				vector<distributed_vector>& velocities,
				const dealii::IndexSet& locally_owned_dofs,
				bool inInitializationProcedure = false) const;

	double getEquilibriumDistribution(size_t i, const numeric_vector& u,
				const double rho) const;*/

	DistributionFunctions& getFormerF()
	{
		return m_formerF;
	}

	DistributionFunctions& getFormerFEq()
	{
		return m_formerFEq;
	}

	void setFormerF(DistributionFunctions& f_tmp)
	{
		m_formerF = f_tmp;
	}

	void setFormerFEq(DistributionFunctions& f_tmp)
	{
		m_formerFEq = f_tmp;
	}

	void initialize(DistributionFunctions& f)
	{
		m_formerF= f;
		m_formerFEq = f;
	}

	/*double getDt() const {
		return m_dt;
	}*/





};

} /* namespace natrium */

#endif /* MULTISTEPCOLLISIONMODEL_H_ */
