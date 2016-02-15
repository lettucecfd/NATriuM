/*
 * Filter.h
 *
 *  Created on: 11.02.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_SMOOTHING_FILTER_H_
#define LIBRARY_NATRIUM_SMOOTHING_FILTER_H_

namespace natrium {

template <size_t dim>
class Filter {
public:
	Filter(){
	}
	virtual ~Filter(){
	}
	virtual void applyFilter(const dealii::DoFHandler<dim>& dof_handler, distributed_vector& dof_vector) = 0;
};

} // namespace natrium

#endif /* LIBRARY_NATRIUM_SMOOTHING_FILTER_H_ */
