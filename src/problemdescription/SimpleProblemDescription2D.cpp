/**
 * @file SimpleProblemDescription2D.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "SimpleProblemDescription2D.h"

namespace natrium {

SimpleProblemDescription2D::SimpleProblemDescription2D(shared_ptr<Triangulation<2> > triangulation, float_t viscosity):
	ProblemDescription(triangulation, viscosity){
}

SimpleProblemDescription2D::~SimpleProblemDescription2D() {
}

} /* namespace natrium */
