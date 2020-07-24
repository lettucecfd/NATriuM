/**
 * @file D3Q45.h
 * @short D3Q45 Stencil
 * @date 02.06.2020
 * @author dwilde3m, University of Siegen, Germany
 */

#ifndef D3Q45MODEL_H_
#define D3Q45MODEL_H_

#include "Stencil.h"

#include "../utilities/BasicNames.h"

namespace natrium {

/** @short D3Q45 Model
 */
class D3Q45: public Stencil {

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
	D3Q45(double scaling = 1.0);

	/// destructor
	virtual ~D3Q45();

	virtual double getSpeedOfSound() const {
		return m_speedOfSound;
	}
	virtual double getSpeedOfSoundSquare() const {
		return m_speedOfSoundSquare;
	}

	virtual size_t getIndexOfOppositeDirection(size_t index) const {
		switch (index){
            case 0 :
                return 0 ;
                break;
            case 1 :
                return 12 ;
                break;
            case 2 :
                return 11 ;
                break;
            case 3 :
                return 10 ;
                break;
            case 4 :
                return 9 ;
                break;
            case 5 :
                return 8 ;
                break;
            case 6 :
                return 7 ;
                break;
            case 7 :
                return 6 ;
                break;
            case 8 :
                return 5 ;
                break;
            case 9 :
                return 4 ;
                break;
            case 10 :
                return 3 ;
                break;
            case 11 :
                return 2 ;
                break;
            case 12 :
                return 1 ;
                break;
            case 13 :
                return 32 ;
                break;
            case 14 :
                return 31 ;
                break;
            case 15 :
                return 30 ;
                break;
            case 16 :
                return 29 ;
                break;
            case 17 :
                return 28 ;
                break;
            case 18 :
                return 27 ;
                break;
            case 19 :
                return 26 ;
                break;
            case 20 :
                return 25 ;
                break;
            case 21 :
                return 24 ;
                break;
            case 22 :
                return 23 ;
                break;
            case 23 :
                return 22 ;
                break;
            case 24 :
                return 21 ;
                break;
            case 25 :
                return 20 ;
                break;
            case 26 :
                return 19 ;
                break;
            case 27 :
                return 18 ;
                break;
            case 28 :
                return 17 ;
                break;
            case 29 :
                return 16 ;
                break;
            case 30 :
                return 15 ;
                break;
            case 31 :
                return 14 ;
                break;
            case 32 :
                return 13 ;
                break;
            case 33 :
                return 44 ;
                break;
            case 34 :
                return 43 ;
                break;
            case 35 :
                return 42 ;
                break;
            case 36 :
                return 41 ;
                break;
            case 37 :
                return 40 ;
                break;
            case 38 :
                return 39 ;
                break;
            case 39 :
                return 38 ;
                break;
            case 40 :
                return 37 ;
                break;
            case 41 :
                return 36 ;
                break;
            case 42 :
                return 35 ;
                break;
            case 43 :
                return 34 ;
                break;
            case 44 :
                return 33 ;
                break;

            default:
				return 0;
		}
	}


	virtual double getMaxParticleVelocityMagnitude() const {
        return sqrt(2)*m_scaling;
	}

	virtual double getScaling() const {
		return m_scaling;
	}

};

} /* namespace natrium */
#endif /* D3Q45MODEL_H_ */
