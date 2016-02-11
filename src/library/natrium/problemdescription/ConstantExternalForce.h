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


template<size_t dim>
class ConstantExternalForce {
private:
	dealii::Tensor<1, dim> m_force;
public:
	ConstantExternalForce(dealii::Tensor<1, dim> force){
		m_force = force;
	}
	virtual ~ConstantExternalForce(){};

	const dealii::Tensor<1, dim>& getForce() const {
		return m_force;
	}

	void scale(double factor) {
		m_force *= factor;
	}
};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_PROBLEMDESCRIPTION_CONSTANTEXTERNALFORCE_H_ */
