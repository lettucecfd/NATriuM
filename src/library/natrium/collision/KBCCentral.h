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

	mutable int counter = 0;
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

	void writeDeviation(double ave, double dev, double ave_entropy, double dev_entropy) const {
		parameterFile << counter << " " << ave << " " << dev << " " << ave_entropy << " " << dev_entropy << endl;
		counter += 1;
	}



private:

	mutable std::ofstream parameterFile;

}
;
} /* namespace natrium */


#endif /* KBCCENTRAL_H_ */
