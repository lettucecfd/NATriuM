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
			const shared_ptr<Stencil> stencil);
	virtual ~MRTStandard();

	virtual void collideAll(DistributionFunctions& f,
			distributed_vector& densities,
			vector<distributed_vector>& velocities,
			bool inInitializationProcedure) const;

	vector<double> m_D;
	vector<double> m_S;

	vector<double> setMRTWeights() const {
		vector<double> D(9);
		D.at(0) = 9;
		D.at(1) = 36;
		D.at(2) = 36;
		D.at(3) = 6;
		D.at(4) = 12;
		D.at(5) = 6;
		D.at(6) = 12;
		D.at(7) = 4;
		D.at(8) = 4;
		return D;
	}

	vector<double> setRelaxationRates() const {
		vector<double> s(9);
		s.at(0) = s.at(3) = s.at(5) = 0;
		s.at(7) = s.at(8) = 1 / getPrefactor();
		s.at(4) = s.at(6) = 8 * (2 - s.at(7)) / (8 - s.at(7));
		s.at(1) = 1.6;
		s.at(2) = 1.8;
		return s;
	}

};

} //namespace natrium

#endif /* MRTSTANDARD_H_ */

