/**
 * @file CouetteFlow2D.h
 * @short Description of a simple Couette Flow (regular channel flow in rectangular domain).
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef COUETTEFLOW2D_H_
#define COUETTEFLOW2D_H_

#include "../problemdescription/SimpleProblemDescription2D.h"

namespace natrium {

/** @short Description of a simple Couette Flow (regular channel flow in rectangular domain).
 */
class CouetteFlow2D: SimpleProblemDescription2D {
public:

	/// constructor
	CouetteFlow2D();

	/// destructor
	virtual ~CouetteFlow2D();
};

} /* namespace natrium */
#endif /* COUETTEFLOW2D_H_ */
