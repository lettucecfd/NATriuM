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

#include "CollisionModel.h"

#include "../utilities/BasicNames.h"

namespace natrium {

/** @short Description of the BGK model for the transformed particle distributions,
 *        as described in Global data which is used by Min and Lee (2011): A spectral-element discontinuous
 *        Galerkin lattice Boltzmann method for nearly incompressible flows, JCP 230 pp. 245-259.
 */
class BGKTransformed: public CollisionModel {

private:

	/// prefactor of the collision (- 1/(tau + 0.5))
	double m_prefactor;

public:

	/**
	 * @short constructor
	 * @param[in] relaxationParameter relaxation parameter tau
	 */
	BGKTransformed(double relaxationParameter, boost::shared_ptr<BoltzmannModel> boltzmannModel);


	/// destructor
	virtual ~BGKTransformed();

	/**
	 * @short function for collision
	 * @param[in/out] distributions the particle distribution functions
	 */
	virtual void collide(vector<double>& distributions) const;


	/**
	 * @short virtual function for collision
	 * @param[in] doF the doF index for which collision is done
	 * @param[in] feq the vector of local equilibrium distributions
	 * @param[in] f the vector of global distribution functions
	 */
	virtual void collide(size_t doF, const vector<double>& feq, vector<distributed_vector>& f);

};

} /* namespace natrium */
#endif /* BGKTRANSFORMED_H_ */
