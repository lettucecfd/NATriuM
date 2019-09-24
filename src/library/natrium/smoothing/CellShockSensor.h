/*
 * CellShockSensor.h
 *
 *  Created on: 16.09.2019
 *      Author: natrium
 */

#ifndef LIBRARY_NATRIUM_SMOOTHING_CELLSHOCKSENSOR_H_
#define LIBRARY_NATRIUM_SMOOTHING_CELLSHOCKSENSOR_H_

#include "../utilities/BasicNames.h"
#include "Filter.h"


namespace natrium {

template<size_t dim>
class CellShockSensor: public Filter<dim> {
private:
	const dealii::FiniteElement<dim>& m_sourceFE;
public:
	CellShockSensor(const dealii::FiniteElement<dim>& fe);
	virtual ~CellShockSensor();
	virtual void applyFilter(const dealii::DoFHandler<dim>& dof_handler,
			distributed_vector& dof_vector);
};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_SMOOTHING_CELLSHOCKSENSOR_H_ */
