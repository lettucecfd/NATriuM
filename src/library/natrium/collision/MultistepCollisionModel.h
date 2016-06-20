/*
 * MultistepCollisionModel.h
 *
 *  Created on: 20.06.2016
 *      Author: dominik
 */

#ifndef MULTISTEPCOLLISIONMODEL_H_
#define MULTISTEPCOLLISIONMODEL_H_

#include <natrium/collision/CollisionModel.h>

namespace natrium {

class MultistepCollisionModel: public CollisionModel {
public:
	MultistepCollisionModel();
	virtual ~MultistepCollisionModel();
};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_COLLISION_MULTISTEPCOLLISIONMODEL_H_ */
