/**
 * @file PhysicalProperties.cpp
 * @short
 * @date 06.06.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "PhysicalProperties.h"

namespace natrium {

template<> PhysicalProperties<2>::PhysicalProperties() {
}
template<> PhysicalProperties<3>::PhysicalProperties() {
}
template<> PhysicalProperties<2>::~PhysicalProperties() {
}
template<> PhysicalProperties<3>::~PhysicalProperties() {
}

template<> double PhysicalProperties<2>::kineticEnergy(
		const vector<distributed_vector>& u, const distributed_vector& rho) {
	const size_t n_dofs = u.at(0).size();
	assert (n_dofs == rho.size());
	assert (n_dofs == u.at(1).size());
	double kE = 0.0;
	for (size_t i = 0; i < n_dofs; i++){
		kE += 0.5 * rho(i) * (u.at(0)(i)*u.at(0)(i) + u.at(1)(i)*u.at(1)(i));
	}
	return kE*0.5/n_dofs;
}
template<> double PhysicalProperties<3>::kineticEnergy(
		const vector<distributed_vector>& u, const distributed_vector& rho) {
	const size_t n_dofs = u.at(0).size();
	assert (n_dofs == rho.size());
	assert (n_dofs == u.at(1).size());
	assert (n_dofs == u.at(2).size());
	double kE = 0.0;
	for (size_t i = 0; i < n_dofs; i++){
		kE += rho(i) * (u.at(0)(i)*u.at(0)(i) + u.at(1)(i)*u.at(1)(i) + u.at(2)(i)*u.at(2)(i));
	}
	return kE*0.5/n_dofs;
}

} /* namespace natrium */
