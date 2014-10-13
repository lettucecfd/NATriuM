/**
 * @file CFDSolverUtilities.h
 * @short General utility functions for the CFD solver
 * @date 04.08.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef CFDSOLVERUTILITIES_H_
#define CFDSOLVERUTILITIES_H_

#include "deal.II/grid/tria.h"
#include "deal.II/base/quadrature_lib.h"

#include "../utilities/BasicNames.h"

namespace natrium {

namespace CFDSolverUtilities {

/**
 * @short
 */
template<size_t dim>
double getMinimumDoFDistanceGLL(const dealii::Triangulation<dim>& tria,
		const size_t orderOfFiniteElement);


/**
 * @short stolen from Deal.II's step 49 tutorial
 */
template<int dim>
void mesh_info(const dealii::Triangulation<dim> &tria,
               const std::string        &filename);

} /* CFDSolverUtilities */
} /* namespace natrium */

#endif /* CFDSOLVERUTILITIES_H_ */
