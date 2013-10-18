/**
 * @file CouetteFlow2D.h
 * @short Description of a simple Couette Flow (regular channel flow in rectangular domain).
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef COUETTEFLOW2D_H_
#define COUETTEFLOW2D_H_

#include "boost/shared_ptr.hpp"
#include "boost/make_shared.hpp"

#include "deal.II/grid/tria.h"

#include "problemdescription/SimpleProblemDescription2D.h"
#include "utilities/BasicNames.h"

using boost::make_shared;
using boost::shared_ptr;
using dealii::Triangulation;


namespace natrium {

/** @short Description of a simple Couette Flow (regular channel flow in rectangular domain).
 */
class CouetteFlow2D: SimpleProblemDescription2D {
public:

	/// constructor
	CouetteFlow2D(float_t viscosity);

	/// destructor
	virtual ~CouetteFlow2D();

private:

	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	shared_ptr<Triangulation<2> > make_grid(){

		//Creation of the principal domain
		shared_ptr<Triangulation<2> > rectMain = make_shared<Triangulation<2> >();

		return rectMain;
	}

};

} /* namespace natrium */
#endif /* COUETTEFLOW2D_H_ */
