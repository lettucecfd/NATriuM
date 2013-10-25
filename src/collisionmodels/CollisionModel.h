/**
 * @file CollisionModel.h
 * @short Abstract class for the description of collision schemes.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef COLLISIONMODEL_H_
#define COLLISIONMODEL_H_

#include "boost/shared_ptr.hpp"

#include "../boltzmannmodels/BoltzmannModel.h"

#include "../utilities/BasicNames.h"


namespace natrium {

/** @short Abstract class for the description of collision schemes.
 */
class CollisionModel {

protected:

	/// relaxation parameter
	double m_relaxationParameter;

	/// Boltzmann model (e.g. D2Q9Incompressible)
	const boost::shared_ptr<BoltzmannModel> m_boltzmannModel;

	/// D (dimension)
	size_t m_d;

	/// Q (number of directions)
	double m_q;

public:

	/**
	 * @short constructor
	 * @param[in] relaxationParameter relaxation parameter tau
	 */
	CollisionModel(double relaxationParameter, boost::shared_ptr<BoltzmannModel> boltzmannModel);

	/// destructor
	virtual ~CollisionModel();

	/**
	 * @short virtual function for collision
	 * @param[in/out] distributions the particle distribution functions
	 */
	virtual void collide(vector<double>& distributions) const = 0;

	/// get relaxation parameter
	double getRelaxationParameter() const {
		return m_relaxationParameter;
	}

};

} /* namespace natrium */
#endif /* COLLISIONMODEL_H_ */
