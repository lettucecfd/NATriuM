/**
 * @file SimpleProblemDescription2D.h
 * @short Description of simple 2D test problems, using boundary IDs and easy-to-use boundary functions.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef SIMPLEPROBLEMDESCRIPTION2D_H_
#define SIMPLEPROBLEMDESCRIPTION2D_H_

#include "ProblemDescription.h"
#include "deal.II/grid/tria.h"

namespace natrium {

/** @short Description of simple 2D test problems, using boundary IDs and easy-to-use boundary functions.
 */
class SimpleProblemDescription2D: public ProblemDescription<2> {
private:

	/// computational grid
	Triangulation<2> m_triangulation;

	/// boundary description

	/// initial velocities

	/// viscosity
	float_t m_viscosity;

public:

	/// constructor
	SimpleProblemDescription2D();

	/// destructor
	virtual ~SimpleProblemDescription2D();
};

} /* namespace natrium */
#endif /* SIMPLEPROBLEMDESCRIPTION2D_H_ */
