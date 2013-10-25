/**
 * @file SimpleProblemDescription2D.h
 * @short Description of simple 2D test problems, using boundary IDs and easy-to-use boundary functions.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef SIMPLEPROBLEMDESCRIPTION2D_H_
#define SIMPLEPROBLEMDESCRIPTION2D_H_

#include "deal.II/grid/tria.h"

#include "../problemdescription/ProblemDescription.h"
#include "../utilities/BasicNames.h"

using dealii::Triangulation;

namespace natrium {

/** @short Description of simple 2D test problems, using boundary IDs and easy-to-use boundary functions.
 */
class SimpleProblemDescription2D: public ProblemDescription<2> {

public:

	/// constructor
	SimpleProblemDescription2D(shared_ptr<Triangulation<2> > triangulation, double viscosity);

	/// destructor
	virtual ~SimpleProblemDescription2D();
};

} /* namespace natrium */
#endif /* SIMPLEPROBLEMDESCRIPTION2D_H_ */
