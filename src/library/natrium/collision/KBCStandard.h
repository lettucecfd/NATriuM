/*
 * KBCStandard.h
 *
 *  Created on: 17.11.2015
 *      Author: Dominik Wilde
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

	double getEquilibriumDistribution(size_t i, const numeric_vector& u,
			const double rho) const;

	/**
	 * @short function for collision
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
	 * @short optimized version of collideAll for D3Q15 stencil
	 */
	void collideAllD3Q15(DistributionFunctions& f,
			distributed_vector& densities,
			vector<distributed_vector>& velocities,
			const dealii::IndexSet& locally_owned_dofs,
			bool inInitializationProcedure) const;

	/**
	 * @short Stores the values of the stabilizer of the KBC function and averages it for every time step
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

	/**
	 * @short writes the averaged data of the stabilizer into a parameter file
	 */

	mutable int counter;


	void writeDeviation(double ave, double dev, double ave_entropy, double dev_entropy) const {
		parameterFile << counter << " " << ave << " " << dev << " " << ave_entropy << " " << dev_entropy << endl;
		counter += 1;
	}




private:

	mutable std::ofstream parameterFile;

}
;
} /* namespace natrium */

#endif /* KBCSTANDARD_H_ */
