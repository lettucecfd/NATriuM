/**
 * @file RungeKutta5LowStorage.h
 * @short Fifth-order Runge-Kutta time integration scheme with low storage consumption.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef RUNGEKUTTA5LOWSTORAGE_H_
#define RUNGEKUTTA5LOWSTORAGE_H_

#include "TimeIntegrator.h"

namespace natrium {

/** @short Implementation of the fifth-order Runge-Kutta time integration scheme with low storage consumption.
 *  @note  The scheme is described in Min and Lee (2011): A spectral-element discontinuous Galerkin lattice
 *         Boltzmann method for nearly incompressible flows, JCP 230 pp. 245-259.
 */
class RungeKutta5LowStorage: public TimeIntegrator {
private:

	/// functions to initialize the RK coefficients
	vector<double> makeA();
	vector<double> makeB();
	vector<double> makeC();

	/// coefficients of the RK scheme
	/// Source: http://www.ece.uvic.ca/~bctill/papers/numacoust/Carpenter_Kennedy_1994.pdf ("solution 3")
	const vector<double> m_a;
	const vector<double> m_b;
	const vector<double> m_c;

	// The scheme has storage 2N (f and df). this vector stores df
	distributed_vector m_df;

	// The storage of F(f) = systemMatrix * f requires an additional storage array
	distributed_vector m_Af;

public:


	/**
	 * @short Constructor
	 * @param timeStepSize The initial time step size
	 * @param problemSize the number of degrees of freedom
	 */
	RungeKutta5LowStorage(double timeStepSize, size_t problemSize);

	/// destructor
	virtual ~RungeKutta5LowStorage();

	/**
	 * @short make one time integration using the RK5 scheme
	 * 			dU_j = A_j * dU_{j-1} + timeStep * F(U_j)
	 * 			U_j  = U_{j-1} + B_j* dU_j
	 */
	virtual void step(distributed_vector& vector,
			const distributed_sparse_matrix& systemMatrix);

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
