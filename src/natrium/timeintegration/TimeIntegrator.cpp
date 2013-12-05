/**
 * @file TimeIntegrator.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "TimeIntegrator.h"

namespace natrium {


vector<vector<double> > natrium::TimeIntegrator::makeA() {
	vector<vector<double> > A(1);
	return A;
}

vector<vector<double> > natrium::TimeIntegrator::makeB() {
	vector<vector<double> > A(1);
		return A;
}

vector<double> natrium::TimeIntegrator::makeC() {
	vector<double> A(1);
		return A;
}

natrium::TimeIntegrator::TimeIntegrator(double timeStepSize):
	m_a(makeA()), m_b(makeB()), m_c(makeC()), m_timeStepSize(timeStepSize){
}

} /* namespace natrium */
