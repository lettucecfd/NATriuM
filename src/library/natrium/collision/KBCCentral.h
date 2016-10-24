/*
 * KBCCentral.h
 *
 *  Created on: 29.03.2016
 *      Author: Dominik Wilde
 */

#ifndef KBCCENTRAL_H_
#define KBCCENTRAL_H_

#include "MRTStandard.h"

#include "../utilities/BasicNames.h"

#include <cassert>

#include "../utilities/Math.h"

namespace natrium {

class KBCCentral: public MRT {
public:
	KBCCentral(double relaxationParameter, double dt,
			const boost::shared_ptr<Stencil> stencil);
	virtual ~KBCCentral();
	virtual void collideAll(DistributionFunctions& f,
			distributed_vector& densities,
			vector<distributed_vector>& velocities,
			const dealii::IndexSet& locally_owned_dofs,
			bool inInitializationProcedure = false) const;

	/**
	 * @short optimized version of collideAll for D2Q9 stencil
	 */
	void collideAllD2Q9(DistributionFunctions& f, distributed_vector& densities,
			vector<distributed_vector>& velocities,
			const dealii::IndexSet& locally_owned_dofs,
			bool inInitializationProcedure) const;

	/**
	 * @short evaluates the stabilizer of the KBC function and writes it to a parameter file
	 */


	struct stabilizer {
		stabilizer(int size_dofs) {
			size = size_dofs;
			value.resize(size);
		}
		int size;
		vector<double> value;
		double average = 0;

		double getDeviation() {
			double deviation = 0;
			average = getAverage();
			for (int r = 0; r < size; r++) {

				deviation += pow(value.at(r) - average, 2);
			}

			deviation /= size;
			deviation = pow(deviation, 0.5);
			return deviation;

		}

		double getAverage() {
			average = 0;
			for (int r = 0; r < size; r++) {
				average += value.at(r);
			}
			average /= size;
			return average;

		}
		;

	};

	//mutable int counter = 0;

//	void writeDeviation(double ave, double dev, double ave_entropy, double dev_entropy) const {
//		parameterFile << counter << " " << ave << " " << dev << " " << ave_entropy << " " << dev_entropy << endl;
//		counter += 1;
//	}

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

		// for D2Q9 only s(1) and s(2) should be adjusted to the specific needs
		vector<double> setRelaxationRates() const {
			vector<double> s(9);
			s.at(0) = s.at(3) = s.at(5) = 0.0;
			s.at(7) = s.at(8) = - getPrefactor();
			s.at(4) = s.at(6) = 8.0 * (2.0 - s.at(7)) / (8.0 - s.at(7));

			s.at(1) = 1.6;
			s.at(2) = 1.8;
			return s;




}

}
;/* namespace natrium */
}

#endif /* KBCCENTRAL_H_ */
