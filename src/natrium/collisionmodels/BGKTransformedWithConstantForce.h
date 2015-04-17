/*

 * BGKTransformedWithConstantForce.h
 *
 *  Created on: Jul 30, 2014
 *      Author: bajat
 */

#ifndef BGKTRANSFORMEDWITHCONSTANTFORCE_H_
#define BGKTRANSFORMEDWITHCONSTANTFORCE_H_

#include "CollisionModel.h"
#include "BGKTransformed.h"
#include "../utilities/BasicNames.h"
#include "../solver/DistributionFunctions.h"
#include "../boltzmannmodels/BoltzmannModel.h"

namespace natrium {

  /** @short Description of the BGK model for the transformed particle distributions with a constant body force acting as the driving force for the flow,
   *
   */
class BGKTransformedWithConstantForce: public BGKTransformed {

private:
  double m_constantBodyForce;
  double m_timeStep;
  distributed_vector m_density;

public:
  BGKTransformedWithConstantForce(double relaxationParameter,
      boost::shared_ptr<BoltzmannModel> boltzmannModel,double timeStep, double constantBodyForce) ;



 void collideSingleDoF(size_t doF, const vector<double>& feq,
      DistributionFunctions& f, distributed_vector& density) const ;
};

}


#endif /* BGKTRANSFORMEDWITHCONSTANTFORCE_H_ */
