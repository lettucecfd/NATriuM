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

	/** coefficients copied from Min's code **/
	/// functions to initialize the RK coefficients
	vector<double> makeA(){
		vector<double> A;
		A += /*0, -0.4178904745, -1.192151694643, -1.697784692471, -1.514183444257;*/
	       0.0,  -567301805773.0/1357537059087.0,  -2404267990393.0/2016746695238.0,  -3550918686646.0/2091501179385.0,  -1275806237668.0/842570457699.0;
		return A;
	}
	vector<double> makeB(){
		vector<double> B;
		B += /*0.1496590219993, 0.3792103129999, 0.8229550293869, 0.6994504559488, 0.1530572479681;*/
			     1432997174477.0/9575080441755.0,   5161836677717.0/13612068292357.0,  1720146321549.0/2090206949498.0,  3134564353537.0/4481467310338.0, 2277821191437.0/14882151754819.0;
		return B;
	}
	vector<double> makeC(){
		vector<double> C;
		C += /*0, 0.1496590219993, 0.3704009573644, 0.6222557631345, 0.9582821306748;*/
				0.0,  1432997174477.0/9575080441755.0,  2526269341429.0/6820363962896.0, 2006345519317.0/3224310063776.0, 2802321613138.0/2924317926251.0,  1.;
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

	/**
	 * @short Constructor
	 * @param timeStepSize The initial time step size
	 * @param problemSize the number of degrees of freedom
	 * @param numberOfBlocks the number of blocks of the block vectors
	 */
	RungeKutta5LowStorage(double timeStepSize, size_t problemSize, size_t numberOfBlocks);

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
