/**
 * @file Stencil.cpp
 * @short Abstract class for the description of a boltzmann model
 * @date 04.06.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "Stencil.h"


namespace natrium {


Stencil::Stencil(size_t d, size_t q,
		 const vector<numeric_vector>& directions,
		 const vector<double>& weights, StencilType stencilType) :
		m_d(d), m_q(q), m_directions(directions), m_weights(weights), m_stencilType(stencilType) {

	// Test if the data is consistent
	assert(d > (size_t) 1);
	assert(d < (size_t) 4);
	assert(q > (size_t) 0);
	assert(directions.size() - q == 0);
	assert(weights.size() - q == 0);
	//for (size_t i = 0; i < q; i++){
	//	assert(directions.at(i).size() - d == 0);
	//}
}

Stencil::~Stencil() {
}


} /* namespace natrium */
