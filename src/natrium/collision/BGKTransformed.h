/**
 * @file BGKTransformed.h
 * @short Description of the BGK model for the transformed particle distributions,
 *        as described in Global data which is used by Min and Lee (2011): A spectral-elemennt discontinuous
 *        Galerkin lattice Boltzmann method for nearly incompressible flows, JCP 230 pp. 245-259.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef BGKTRANSFORMED_H_
#define BGKTRANSFORMED_H_

#include "Collision.h"

#include "../utilities/BasicNames.h"
#include "../solver/DistributionFunctions.h"

namespace natrium {

/** @short Description of the BGK model for the transformed particle distributions,
 *        as described in Global data which is used by Min and Lee (2011): A spectral-element discontinuous
 *        Galerkin lattice Boltzmann method for nearly incompressible flows, JCP 230 pp. 245-259.
 */
class BGKTransformed {
private:

	// number of discrete particle velocities
	size_t m_Q;

	/// relaxation parameter
	double m_relaxationParameter;

	/// prefactor of the collision (- 1/(tau + 0.5))
	double m_prefactor;

public:

	/**
	 * @short constructor
	 * @param[in] relaxationParameter relaxation parameter tau
	 */
	BGKTransformed(size_t Q, double relaxationParameter);

	/// destructor
	virtual ~BGKTransformed();

	/* @short virtual function for collision
	* @param[in] doF the doF index for which collision is done
	* @param[in] feq the vector of local equilibrium distributions
	* @param[in] f the vector of global distribution functions
	*/
	void collideSingleDoF(size_t doF, const vector<double>& feq,
			DistributionFunctions& f) const {
		for (size_t j = 0; j < m_Q; j++) {
			f.at(j)(doF) += m_prefactor * (f.at(j)(doF) - feq.at(j));
		}

	}
	/* collide */

	/**
	 * @short calculate relaxation parameter
	 */
	static double calculateRelaxationParameter(double viscosity,
			double timeStepSize,
			const BoltzmannModel& boltzmannModel) {
		assert(viscosity > 0.0);
		assert(timeStepSize > 0.0);
		return (viscosity)
				/ (timeStepSize * boltzmannModel.getSpeedOfSoundSquare());
	}

	double getPrefactor() const {
		return m_prefactor;
	}

	size_t getQ() const {
		return m_Q;
	}

	double getRelaxationParameter() const {
		return m_relaxationParameter;
	}
};

} /* namespace natrium */
#endif /* BGKTRANSFORMED_H_ */
