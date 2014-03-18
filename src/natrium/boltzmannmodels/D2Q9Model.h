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
	vector<numeric_vector> makeDirections(double scaling);

	/** @short function to create the vector of weights
	 *  @return the vector of weights
	 */
	vector<double> makeWeights();

protected:

	/// speed of sound
	const double m_speedOfSound;

	/// (speed of sound)^2
	const double m_speedOfSoundSquare;

public:

	/// D
	static const size_t D;

	/// Q
	static const size_t Q;

	/// constructor
	D2Q9Model(double scaling);

	/// destructor
	virtual ~D2Q9Model();

	virtual double getSpeedOfSound() const {
		return m_speedOfSound;
	}
	virtual double getSpeedOfSoundSquare() const {
		return m_speedOfSoundSquare;
	}

	virtual size_t getIndexOfOppositeDirection(size_t index)const {
		switch (index){
			case 0:
				return 0;
				break;
			case 1:
				return 3;
				break;
			case 2:
				return 4;
				break;
			case 3:
				return 1;
				break;
			case 4:
				return 2;
				break;
			case 5:
				return 7;
				break;
			case 6:
				return 8;
				break;
			case 7:
				return 5;
				break;
			case 8:
				return 6;
				break;
			default:
				return 0;
		}
	}

};

} /* namespace natrium */
#endif /* D2Q9MODEL_H_ */
