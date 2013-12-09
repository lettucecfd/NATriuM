/**
 * @file TimeIntegrator.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "TimeIntegrator.h"

namespace natrium {


TimeIntegrator::TimeIntegrator(double timeStepSize):
	m_timeStepSize(timeStepSize){
}

/// destructor
TimeIntegrator::~TimeIntegrator(){
}

} /* namespace natrium */
