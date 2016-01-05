/*
 * ConstantExternalForce.h
 *
 *  Created on: 04.01.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_PROBLEMDESCRIPTION_CONSTANTEXTERNALFORCE_H_
#define LIBRARY_NATRIUM_PROBLEMDESCRIPTION_CONSTANTEXTERNALFORCE_H_

#include "deal.II/base/tensor.h"

namespace natrium {

enum ForceType {
	NO_FORCING, 	// No external force
	SHIFTING_VELOCITY,  // Shifting velocity method
	EXACT_DIFFERENCE,   // Exact difference method by Kuppershtokh
	GUO
};

template<size_t dim>
class ConstantExternalForce {
private:
	dealii::Tensor<1, dim> m_force;
	ForceType m_forceType;
public:
	ConstantExternalForce(dealii::Tensor<1, dim> force,
	ForceType force_type){
		m_force = force;
		m_forceType = force_type;
	}
	virtual ~ConstantExternalForce();

	const dealii::Tensor<1, dim>& getForce() const {
		return m_force;
	}

	ForceType getForceType() const {
		return m_forceType;
	}
};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_PROBLEMDESCRIPTION_CONSTANTEXTERNALFORCE_H_ */
