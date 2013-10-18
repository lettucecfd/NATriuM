/**
 * @file ProblemDescription.h
 * @short Abstract class for the description of a CFD problem. The description includes the computational mesh,
 *        boundary description, viscosity and initial values.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef PROBLEMDESCRIPTION_H_
#define PROBLEMDESCRIPTION_H_

#include "boost/shared_ptr.hpp"

#include "deal.II/grid/tria.h"

#include "../utilities/BasicNames.h"

using dealii::Triangulation;
using boost::shared_ptr;

namespace natrium {

/** @short Abstract class for the description of a CFD problem. The description includes the computational mesh,
 *         boundary description, viscosity and initial values.
 *  @tparam dim The dimension of the flow (2 or 3).
 */
template<int dim> class ProblemDescription {
private:

	/// computational grid
	shared_ptr<Triangulation<2> > m_triangulation;

	/// boundary description

	/// initial velocities

	/// viscosity
	float_t m_viscosity;

public:

	/// constructor
	ProblemDescription(shared_ptr<Triangulation<2> > triangulation, float_t viscosity);

	///  destructor
	virtual ~ProblemDescription(){};
};

} /* namespace natrium */


#endif /* PROBLEMDESCRIPTION_H_ */
