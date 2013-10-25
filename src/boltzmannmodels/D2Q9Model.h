/**
 * @file D2Q9Model.h
 * @short D2Q9 Boltzmann Model
 * @date 02.06.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef D2Q9MODEL_H_
#define D2Q9MODEL_H_

#include "BoltzmannModel.h"

#include "../utilities/BasicNames.h"

namespace natrium {

/** @short D2Q9 Model
 */
class D2Q9Model: public BoltzmannModel {

private:

	/** @short function to create the vector of directions
	 *  @return the vector of directions
	 */
	vector<numeric_vector> makeDirections();

	/** @short function to create the vector of weights
	 *  @return the vector of weights
	 */
	vector<double> makeWeights();

public:

	/// D
	static const size_t D;

	/// Q
	static const size_t Q;

	/// speed of sound
	static const double speedOfSound;

	/// (speed of sound)^2
	static const double speedOfSoundSquare;

	/// constructor
	D2Q9Model();

	/// destructor
	virtual ~D2Q9Model();
};

} /* namespace natrium */
#endif /* D2Q9MODEL_H_ */
