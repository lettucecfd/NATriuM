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

#include "../stencils/Stencil.h"

#include "../utilities/BasicNames.h"
#include "../utilities/NATriuMException.h"

namespace natrium {

namespace CFDSolverUtilities {

class CFDSolverUtilitiesException: public NATriuMException {
private:
	std::string message;
public:
	CFDSolverUtilitiesException(const char *msg) :
		NATriuMException(msg), message(msg) {
	}
	CFDSolverUtilitiesException(const string& msg) :
		NATriuMException(msg), message(msg) {
	}
	~CFDSolverUtilitiesException() throw () {
	}
	const char *what() const throw () {
		return this->message.c_str();
	}
};

/**
 * @short
 */
template<size_t dim>
double getMinimumDoFDistanceGLL(const Mesh<dim>& tria,
		const size_t orderOfFiniteElement);

template<size_t dim>
double getMinimumVertexDistance(const Mesh<dim>& tria);

template<size_t dim>
double calculateTimestep(const Mesh<dim>& tria,
		const size_t orderOfFiniteElement, const Stencil& stencil,
		double cFL = 0.4);

/**
 * @short stolen from Deal.II's step 49 tutorial
 */
template<int dim>
void mesh_info(const Mesh<dim> &tria,
		const std::string &filename);

} /* CFDSolverUtilities */
} /* namespace natrium */

#endif /* CFDSOLVERUTILITIES_H_ */
