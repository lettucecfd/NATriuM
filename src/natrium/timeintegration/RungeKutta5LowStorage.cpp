/**
 * @file RungeKutta5LowStorage.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "RungeKutta5LowStorage.h"

#include <cassert>

// enable + operator for filling vectors
#include "boost/assign/std/vector.hpp"
using namespace boost::assign;

namespace natrium {

vector<double> RungeKutta5LowStorage::makeA() {
	vector<double> A;
	A += 0, -0.4178904745, -1.192151694643, -1.697784692471, -1.514183444257;
	return A;
}

vector<double> RungeKutta5LowStorage::makeB() {
	vector<double> B;
	B += 0.1496590219993, 0.3792103129999, 0.8229550293869, 0.6994504559488, 0.1530572479681;
	return B;
}

vector<double> natrium::RungeKutta5LowStorage::makeC() {
	vector<double> C;
	C += 0, 0.1496590219993, 0.3704009573644, 0.6222557631345, 0.9582821306748;
	return C;
}

RungeKutta5LowStorage::RungeKutta5LowStorage(double timeStepSize,
		size_t problemSize): TimeIntegrator(timeStepSize), m_a(makeA()), m_b(makeB()), m_c(makeC()), m_df(
				problemSize), m_Af(problemSize) {

}

RungeKutta5LowStorage::~RungeKutta5LowStorage() {
}

void RungeKutta5LowStorage::step(distributed_vector& f,
		const distributed_sparse_matrix& systemMatrix) {
	// Test all dimensions and change, if necessary
	assert(systemMatrix.n() == systemMatrix.m());
	assert(f.size() == systemMatrix.n());
	if (m_Af.size() != f.size()) {
		m_Af.reinit(f.size());
	}
	if (m_df.size() != f.size()) {
		m_df.reinit(f.size());
	}
	// According to Carpenter and Kennedy (1994) - p. 3
	// df = a*df + h*Af
	// f = f + B* df
	// make first step manually
	systemMatrix.vmult(m_Af,f);
	m_Af *= getTimeStepSize();
	m_df = m_Af;
	m_df *= m_b.at(0);
	f += m_df;
	m_df /= m_b.at(0);
	for (size_t i = 1; i < 5; i++){
		systemMatrix.vmult(m_Af,f);
		m_Af *= getTimeStepSize();
		m_df *= m_a.at(i);
		m_df += m_Af;
		m_df *= m_b.at(i);
		f += m_df;
		m_df /= m_b.at(i);
	}

}

} /* namespace natrium */
