/*
 * MRTStandard.h
 *
 *  Created on: 25.07.2015
 *      Author: dominik
 */

#ifndef MRTSTANDARD_H_
#define MRTSTANDARD_H_

#include "MRT.h"

#include "../utilities/BasicNames.h"

#include <cassert>

#include "../utilities/Math.h"

namespace natrium {

class MRTStandard: public MRT {
public:
	MRTStandard(double relaxationParameter, double dt,
			const boost::shared_ptr<Stencil> stencil);
	virtual ~MRTStandard();

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
			const dealii::IndexSet& locally_owned_dofs,
			bool inInitializationProcedure = false) const;

	/**
	 * @short optimized version of collideAll for D2Q9 stencil
	 */
	void collideAllD2Q9(DistributionFunctions& f,
			distributed_vector& densities, vector<distributed_vector>& velocities,
			const dealii::IndexSet& locally_owned_dofs,
			bool inInitializationProcedure) const;


	vector<double> m_D;
	vector<double> m_S;

	vector<double> setMRTWeights() const {
		vector<double> D(9);
		D.at(0) = 9.0;
		D.at(1) = 36.0;
		D.at(2) = 36.0;
		D.at(3) = 6.0;
		D.at(4) = 12.0;
		D.at(5) = 6.0;
		D.at(6) = 12.0;
		D.at(7) = 4.0;
		D.at(8) = 4.0;
		return D;
	}

	vector<double> setRelaxationRates() const {
		vector<double> s(9);
		s.at(0) = s.at(3) = s.at(5) = 0.0;
		s.at(7) = s.at(8) = - 1.0 / getPrefactor();
		s.at(4) = s.at(6) = 8.0 * (2.0 - s.at(7)) / (8.0 - s.at(7));
		s.at(1) = - 1.0 / getPrefactor();
		s.at(2) = - 1.0 / getPrefactor();
		return s;
	}

};

} //namespace natrium

#endif /* MRTSTANDARD_H_ */

