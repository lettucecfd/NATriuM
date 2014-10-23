/**
 * @file BGKTransformed.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "BGKTransformed.h"

namespace natrium {

/// constructor
BGKTransformed::BGKTransformed(size_t Q, double relaxationParameter) :
	 m_Q(Q), m_relaxationParameter(relaxationParameter), m_prefactor(
				-1. / (relaxationParameter + 0.5)) {


} // constructor

// destructor
BGKTransformed::~BGKTransformed() {
} // destructor


}
/* namespace natrium */

