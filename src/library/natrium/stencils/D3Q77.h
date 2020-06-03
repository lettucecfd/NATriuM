/**
 * @file D3Q77.h
 * @short D3Q77 Stencil
 * @date 02.06.2020
 * @author dwilde3m, University of Siegen, Germany
 */

#ifndef D3Q77MODEL_H_
#define D3Q77MODEL_H_

#include "Stencil.h"

#include "../utilities/BasicNames.h"

namespace natrium {

/** @short D3Q77 Model
 */
class D3Q77: public Stencil {

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
	D3Q77(double scaling = 1.0);

	/// destructor
	virtual ~D3Q77();

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
                return 2 ;
                break;
            case 2 :
                return 1 ;
                break;
            case 3 :
                return 4 ;
                break;
            case 4 :
                return 3 ;
                break;
            case 5 :
                return 6 ;
                break;
            case 6 :
                return 5 ;
                break;
            case 7 :
                return 8 ;
                break;
            case 8 :
                return 7 ;
                break;
            case 9 :
                return 10 ;
                break;
            case 10 :
                return 9 ;
                break;
            case 11 :
                return 12 ;
                break;
            case 12 :
                return 11 ;
                break;
            case 13 :
                return 16 ;
                break;
            case 14 :
                return 15 ;
                break;
            case 15 :
                return 14 ;
                break;
            case 16 :
                return 13 ;
                break;
            case 17 :
                return 20 ;
                break;
            case 18 :
                return 19 ;
                break;
            case 19 :
                return 18 ;
                break;
            case 20 :
                return 17 ;
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
                return 28 ;
                break;
            case 26 :
                return 27 ;
                break;
            case 27 :
                return 26 ;
                break;
            case 28 :
                return 25 ;
                break;
            case 29 :
                return 32 ;
                break;
            case 30 :
                return 31 ;
                break;
            case 31 :
                return 30 ;
                break;
            case 32 :
                return 29 ;
                break;
            case 33 :
                return 36 ;
                break;
            case 34 :
                return 35 ;
                break;
            case 35 :
                return 34 ;
                break;
            case 36 :
                return 33 ;
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
                return 44 ;
                break;
            case 42 :
                return 43 ;
                break;
            case 43 :
                return 42 ;
                break;
            case 44 :
                return 41 ;
                break;
            case 45 :
                return 48 ;
                break;
            case 46 :
                return 47 ;
                break;
            case 47 :
                return 46 ;
                break;
            case 48 :
                return 45 ;
                break;
            case 49 :
                return 52 ;
                break;
            case 50 :
                return 51 ;
                break;
            case 51 :
                return 50 ;
                break;
            case 52 :
                return 49 ;
                break;
            case 53 :
                return 56 ;
                break;
            case 54 :
                return 55 ;
                break;
            case 55 :
                return 54 ;
                break;
            case 56 :
                return 53 ;
                break;
            case 57 :
                return 60 ;
                break;
            case 58 :
                return 59 ;
                break;
            case 59 :
                return 58 ;
                break;
            case 60 :
                return 57 ;
                break;
            case 61 :
                return 68 ;
                break;
            case 62 :
                return 67 ;
                break;
            case 63 :
                return 66 ;
                break;
            case 64 :
                return 65 ;
                break;
            case 65 :
                return 64 ;
                break;
            case 66 :
                return 63 ;
                break;
            case 67 :
                return 62 ;
                break;
            case 68 :
                return 61 ;
                break;
            case 69 :
                return 76 ;
                break;
            case 70 :
                return 75 ;
                break;
            case 71 :
                return 74 ;
                break;
            case 72 :
                return 73 ;
                break;
            case 73 :
                return 72 ;
                break;
            case 74 :
                return 71 ;
                break;
            case 75 :
                return 70 ;
                break;
            case 76 :
                return 69 ;
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
#endif /* D3Q77MODEL_H_ */
