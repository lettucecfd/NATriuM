/*
 * DataProcessor.h
 *
 *  Created on: 05.03.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_DATAPROCESSORS_DATAPROCESSOR_H_
#define LIBRARY_NATRIUM_DATAPROCESSORS_DATAPROCESSOR_H_

#include "../utilities/BasicNames.h"

namespace natrium {

// forward declaration
template <size_t dim>
class CFDSolver;


template <size_t dim>
class DataProcessor{
protected:
	const CFDSolver<3> & m_solver;
public:
	DataProcessor(const CFDSolver<3> & solver):
		m_solver(solver){
	}
	virtual ~DataProcessor(){
	}
	virtual void apply() = 0;
};


} /* namespace natrium */


#endif /* LIBRARY_NATRIUM_DATAPROCESSORS_DATAPROCESSOR_H_ */
