/*

 * BGKTransformedWithConstantForce.cpp
 *
 *  Created on: Jul 30, 2014
 *      Author: bajat
 */
#include "BGKTransformedWithConstantForce.h"
namespace natrium {

  /// constructor
  BGKTransformedWithConstantForce::BGKTransformedWithConstantForce(double relaxationParameter,
      boost::shared_ptr<BoltzmannModel> boltzmannModel , double timeStep, double constantBodyForce):
         BGKTransformed(relaxationParameter, boltzmannModel),
         m_constantBodyForce(constantBodyForce),
         m_timeStep(timeStep)

        {

        } // constructor

  // destructor



  void BGKTransformedWithConstantForce::collideSingleDoF(size_t doF, const vector<double>& feq,
      DistributionFunctions& f, distributed_vector& density) const {
      numeric_vector force(m_d);
      force(0) = m_constantBodyForce;

    for (size_t j = 0; j < m_q; j++) {

      double forceComponent = 3*m_boltzmannModel->getDirection(j)(0)*force(0);
      forceComponent *= forceComponent * m_boltzmannModel->getWeight(j)*density(doF);


      f.at(j)(doF) += m_prefactor * (f.at(j)(doF) - feq.at(j)) + m_timeStep*forceComponent;
    }


}
}


