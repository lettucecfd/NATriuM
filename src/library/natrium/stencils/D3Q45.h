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
	static constexpr size_t D = 3;

	/// Q
	static constexpr size_t Q = 45;

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
            case 1 :
                return 12 ;
            case 2 :
                return 11 ;
            case 3 :
                return 10 ;
            case 4 :
                return 9 ;
            case 5 :
                return 8 ;
            case 6 :
                return 7 ;
            case 7 :
                return 6 ;
            case 8 :
                return 5 ;
            case 9 :
                return 4 ;
            case 10 :
                return 3 ;
            case 11 :
                return 2 ;
            case 12 :
                return 1 ;
            case 13 :
                return 32 ;
            case 14 :
                return 31 ;
            case 15 :
                return 30 ;
            case 16 :
                return 29 ;
            case 17 :
                return 28 ;
            case 18 :
                return 27 ;
            case 19 :
                return 26 ;
            case 20 :
                return 25 ;
            case 21 :
                return 24 ;
            case 22 :
                return 23 ;
            case 23 :
                return 22 ;
            case 24 :
                return 21 ;
            case 25 :
                return 20 ;
            case 26 :
                return 19 ;
            case 27 :
                return 18 ;
            case 28 :
                return 17 ;
            case 29 :
                return 16 ;
            case 30 :
                return 15 ;
            case 31 :
                return 14 ;
            case 32 :
                return 13 ;
            case 33 :
                return 44 ;
            case 34 :
                return 43 ;
            case 35 :
                return 42 ;
            case 36 :
                return 41 ;
            case 37 :
                return 40 ;
            case 38 :
                return 39 ;
            case 39 :
                return 38 ;
            case 40 :
                return 37 ;
            case 41 :
                return 36 ;
            case 42 :
                return 35 ;
            case 43 :
                return 34 ;
            case 44 :
                return 33 ;
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
