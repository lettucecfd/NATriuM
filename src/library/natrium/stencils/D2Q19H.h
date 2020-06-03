/**
 * @file D2Q19H.h
 * @short D2Q19H Stencil
 * @date 02.06.2020
 * @author dwilde3m, University of Siegen, Germany
 */

#ifndef D2Q19HMODEL_H_
#define D2Q19HMODEL_H_

#include "Stencil.h"

#include "../utilities/BasicNames.h"

namespace natrium {

/** @short D2Q19H Model
 */
class D2Q19H: public Stencil {

private:

	/** @short function to create the vector of directions
	 *  @return the vector of directions
	 */
	vector<numeric_vector> makeDirections(double scaling);

	/** @short function to create the vector of weights
	 *  @return the vector of weights
	 */
	vector<double> makeWeights();

	numeric_matrix makeMomentBasis(vector<numeric_vector> e);

protected:

	/// speed of sound
	const double m_speedOfSound;

	/// (speed of sound)^2
	const double m_speedOfSoundSquare;

	/// scaling of the stencil
	const double m_scaling;

public:

	/// D
	static const size_t D;

	/// Q
	static const size_t Q;

	/// constructor
	D2Q19H(double scaling = 1.0);

	/// destructor
	virtual ~D2Q19H();

	virtual double getSpeedOfSound() const {
		return m_speedOfSound;
	}
	virtual double getSpeedOfSoundSquare() const {
		return m_speedOfSoundSquare;
	}

	virtual size_t getIndexOfOppositeDirection(size_t index) const {
		switch (index){
			case 0:
				return 0;
				break;
			case 1:
				return 4;
				break;
			case 2:
				return 3;
				break;
			case 3:
				return 2;
				break;
			case 4:
				return 1;
				break;
			case 5:
				return 8;
				break;
			case 6:
				return 7;
				break;
			case 7:
				return 6;
				break;
			case 8:
				return 5;
				break;
			case 9:
				return 12;
				break;
			case 10:
				return 11;
				break;
			case 11:
				return 10;
				break;
			case 12:
				return 9;
				break;
			case 13:
				return 14;
				break;
			case 14:
				return 13;
				break;
			case 15:
				return 16;
				break;
			case 16:
				return 15;
				break;
			case 17:
				return 18;
				break;
			case 18:
				return 17;
				break;
			/*case 19:
				return 17;
				break;
			case 20:
				return 18;
				break;
			case 21:
				return 23;
				break;
			case 22:
				return 24;
				break;
			case 23:
				return 21;
				break;
			case 24:
				return 22;
				break; */
			default:
				return 0;
		}
	}

    const double m_directionsArray[19][2] {{ 0.0 , 0.0 },{ 1.367469636752619 , 0.775196278121181 },
                                         { 1.367469636752619 , -0.775196278121181 },
                                         { -1.367469636752619 , 0.775196278121181 },
                                         { -1.367469636752619 , -0.775196278121181 },
                                         { 2.6987507639352253 , 1.8663975507141328 },
                                         { 2.6987507639352253 , -1.8663975507141328 },
                                         { -2.6987507639352253 , 1.8663975507141328 },
                                         { -2.6987507639352253 , -1.8663975507141328 },
                                         { 1.105629214668943 , 2.5175897644357486 },
                                         { 1.105629214668943 , -2.5175897644357486 },
                                         { -1.105629214668943 , 2.5175897644357486 },
                                         { -1.105629214668943 , -2.5175897644357486 },
                                         { 2.9213306655318734 , 0.0 },
                                         { -2.9213306655318734 , 0.0 },
                                         { 0.0 , 1.4869982213169028 },
                                         { 0.0 , -1.4869982213169028 },
                                         { 0.0 , 3.8358342053914734 },
                                         { 0.0 , -3.8358342053914734 }};

	virtual double getMaxParticleVelocityMagnitude() const {
        return sqrt(2)*m_scaling;
	}

	virtual double getScaling() const {
		return m_scaling;
	}

};

} /* namespace natrium */
#endif /* D2Q19HMODEL_H_ */
