/**
 * @file ProblemDescription.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "ProblemDescription.h"

namespace natrium {

template<int dim> ProblemDescription<dim>::ProblemDescription(
		shared_ptr<Triangulation<2> > triangulation, float_t viscosity):
			m_triangulation(triangulation),
			m_viscosity(viscosity){
}

} /* namespace natrium */
