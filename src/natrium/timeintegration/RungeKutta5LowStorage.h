/**
 * @file RungeKutta5LowStorage.h
 * @short Fifth-order Runge-Kutta time integration scheme with low storage consumption.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef RUNGEKUTTA5LOWSTORAGE_H_
#define RUNGEKUTTA5LOWSTORAGE_H_

#include "TimeIntegrator.h"
#include "../utilities/BasicNames.h"

namespace natrium {

/** @short Implementation of the fifth-order Runge-Kutta time integration scheme with low storage consumption.
 *  @note  The scheme is described in Min and Lee (2011): A spectral-element discontinuous Galerkin lattice
 *         Boltzmann method for nearly incompressible flows, JCP 230 pp. 245-259.
 */
template <class MATRIX, class VECTOR> class RungeKutta5LowStorage: public TimeIntegrator<MATRIX, VECTOR> {
private:

	/// functions to initialize the RK coefficients
	vector<double> makeA(){
		vector<double> A;
		A += 0, -0.4178904745, -1.192151694643, -1.697784692471, -1.514183444257;
		return A;
	}
	vector<double> makeB(){
		vector<double> B;
		B += 0.1496590219993, 0.3792103129999, 0.8229550293869, 0.6994504559488, 0.1530572479681;
		return B;
	}
	vector<double> makeC(){
		vector<double> C;
		C += 0, 0.1496590219993, 0.3704009573644, 0.6222557631345, 0.9582821306748;
		return C;
	}

	/// coefficients of the RK scheme
	/// Source: http://www.ece.uvic.ca/~bctill/papers/numacoust/Carpenter_Kennedy_1994.pdf ("solution 3")
	const vector<double> m_a;
	const vector<double> m_b;
	const vector<double> m_c;

	// The scheme has storage 2N (f and df). this vector stores df
	VECTOR m_df;

	// The storage of F(f) = systemMatrix * f requires an additional storage array
	VECTOR m_Af;

public:


	/**
	 * @short Constructor
	 * @param timeStepSize The initial time step size
	 * @param problemSize the number of degrees of freedom
	 */
	RungeKutta5LowStorage(double timeStepSize, size_t problemSize);

	/// destructor
	virtual ~RungeKutta5LowStorage(){};

	/**
	 * @short make one time integration using the RK5 scheme
	 * 			dU_j = A_j * dU_{j-1} + timeStep * F(U_j)
	 * 			U_j  = U_{j-1} + B_j* dU_j
	 */
	virtual void step(VECTOR& vector,
			const MATRIX& systemMatrix);

	const vector<double>& getA() const {
		return m_a;
	}

	const vector<double>& getB() const {
		return m_b;
	}

	const vector<double>& getC() const {
		return m_c;
	}
};

} /* namespace natrium */
#endif /* RUNGEKUTTA5LOWSTORAGE_H_ */
