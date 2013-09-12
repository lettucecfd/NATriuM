/**
 * @file CollisionModel.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "CollisionModel.h"

namespace natrium {

// constructor
CollisionModel::CollisionModel(float_t relaxationParameter,
		const boost::shared_ptr<BoltzmannModel> boltzmannModel) :
		m_relaxationParameter(relaxationParameter),
		m_boltzmannModel(boltzmannModel),
		m_d(boltzmannModel->getD()),
		m_q(boltzmannModel->getQ()){

} // constructor

CollisionModel::~CollisionModel() {
}

} /* namespace natrium */
