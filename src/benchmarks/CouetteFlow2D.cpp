/**
 * @file CouetteFlow2D.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "CouetteFlow2D.h"

namespace natrium {


CouetteFlow2D::CouetteFlow2D(float_t viscosity):
	SimpleProblemDescription2D(make_grid(), viscosity){
}

CouetteFlow2D::~CouetteFlow2D() {
}

} /* namespace natrium */
