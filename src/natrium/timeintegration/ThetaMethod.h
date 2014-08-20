/**
 * @file ThetaMethod.h
 * @short Time integration scheme: theta = 0: Explicit Euler, theta = 1: implicit Euler, theta = 0.5: Crank-Nicholson.
 * @date 19.08.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef THETAMETHOD_H_
#define THETAMETHOD_H_

#include "TimeIntegrator.h"
#include "../utilities/BasicNames.h"

namespace natrium {

/** @short Implementation of the Theta method for time integration.
 */
template <class MATRIX, class VECTOR> class ThetaMethod: public TimeIntegrator<MATRIX, VECTOR> {
private:

	double m_theta;

	MATRIX m_tmpMatrix;

	VECTOR m_tmpSystemVector;

public:


	/**
	 * @short Constructor
	 * @param timeStepSize The initial time step size
	 * @param problemSize the number of degrees of freedom
	 * @param theta  theta = 0: Explicit Euler, theta = 1: implicit Euler, theta = 0.5: Crank-Nicholson.
	 */
	ThetaMethod(double timeStepSize, size_t problemSize, double theta);

	/**
	 * @short Constructor
	 * @param timeStepSize The initial time step size
	 * @param problemSize the number of degrees of freedom
	 * @param theta  theta = 0: Explicit Euler, theta = 1: implicit Euler, theta = 0.5: Crank-Nicholson.
	 */
	ThetaMethod(double timeStepSize, size_t problemSize, size_t numberOfBlocks, double theta);

	/// destructor
	virtual ~ThetaMethod(){};

	/**
	 * @short One time step of the theta method
	 */
	 void step(VECTOR& vector,
			const MATRIX& systemMatrix, const VECTOR& systemVector);

};

} /* namespace natrium */
#endif /* RUNGEKUTTA5LOWSTORAGE_H_ */
