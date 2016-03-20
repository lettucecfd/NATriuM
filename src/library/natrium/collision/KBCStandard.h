/*
 * KBCStandard.h
 *
 *  Created on: 17.11.2015
 *      Author: dominik
 */

#ifndef KBCSTANDARD_H_
#define KBCSTANDARD_H_

#include "MRT.h"

#include "../utilities/BasicNames.h"

#include <cassert>

#include "../utilities/Math.h"

namespace natrium {

class KBCStandard: public MRT {
public:
	KBCStandard(double relaxationParameter, double dt,
			const boost::shared_ptr<Stencil> stencil);
	virtual ~KBCStandard();

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
	void collideAllD2Q9(DistributionFunctions& f, distributed_vector& densities,
			vector<distributed_vector>& velocities,
			const dealii::IndexSet& locally_owned_dofs,
			bool inInitializationProcedure) const;

	/**
	 * @short optimized version of collideAll for D2Q9 stencil
	 */
	void collideAllD3Q15(DistributionFunctions& f,
			distributed_vector& densities,
			vector<distributed_vector>& velocities,
			const dealii::IndexSet& locally_owned_dofs,
			bool inInitializationProcedure) const;

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

	mutable int counter;

private:

	mutable std::ofstream parameterFile;

}
;
} /* namespace natrium */

#endif /* KBCSTANDARD_H_ */
