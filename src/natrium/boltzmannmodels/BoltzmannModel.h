/**
 * @file BoltzmannModel.h
 * @short Abstract class for the description of a boltzmann model
 * @date 04.06.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef BOLTZMANNMODEL_H_
#define BOLTZMANNMODEL_H_

#include "assert.h"

#include "../utilities/Math.h"
#include "../utilities/BasicNames.h"

namespace natrium {

/**
 * @short Enum type for the difference stencil
 */
enum StencilType {

	/// D2Q9 stencil
	Stencil_D2Q9

};

/**
 * @short Abstract class for the description of a boltzmann model, e.g. D2Q9
 */
class BoltzmannModel {

private:
	/// dimension
	const size_t m_d;

	/// number of directions
	const size_t m_q;

	/// directions
	const vector<numeric_vector> m_directions;

	/// weights
	const vector<double> m_weights;

	/// DdQq type (id)
	const StencilType m_stencilType;

public:

//////////////////////////////////
// CONSTRUCTION AND DESTRUCTION //
//////////////////////////////////

	/** @short constructor
	 *  @param d dimension
	 *  @param q number of directions
	 *  @param directions the directions of the stencil
	 *  @param weights the weights of the equilibrium distribution
	 *  @param stencilType type of the stencil (e.g. D2Q9)
	 */
	BoltzmannModel(size_t d, size_t q, const vector<numeric_vector>& directions,
			const vector<double>& weights, StencilType stencilType);

	/// destructor
	virtual ~BoltzmannModel();



///////////////////////
// GETTER AND SETTER //
///////////////////////

	/** @short get a reference to the vector of directions
	 *  @return a ublas_vector, which contains the directions of the DdQq-stencil as ublas_vectors
	 */
	const vector<numeric_vector>& getDirections() const {
		return m_directions;
	}

	/** @short get the i-th direction
	 *  @param i index i
	 *  @return a reference to the i-th direction of the DdQq-stencil
	 */
	const numeric_vector& getDirection(size_t i) const {
		assert(i < m_q);
		return m_directions.at(i);
	}

	/** @short get q, the number of directions in the DdQq-stencil
	 *  @return q
	 */
	size_t getQ() const {
		return m_q;
	}

	/** @short get d, the dimension of the DdQq-stencil
	 *  @return d
	 */
	size_t getD() const {
		return m_d;
	}

	/** @short get the weights of the equilibrium distributions
	 *  @return a reference to the vector of weights
	 */
	const vector<double>& getWeights() const {
		return m_weights;
	}

	/** @short get the weight belonging to a certain direction
	 *  @param i index i of the direction (1 <= i <= q)
	 *  @return the i-th weight
	 */
	double getWeight(size_t i) const {
		assert(i < m_q);
		return m_weights.at(i);
	}

	/** @short get stencil type
	 *  @return stencil type, e.g. D2Q9
	 */
	const StencilType getStencilType() const {
		return m_stencilType;
	}


////////////////////////////////////
// CALCULATE MACROSCOPIC ENTITIES //
////////////////////////////////////

	/**
	 * @short calculate macroscopic density
	 * @param[in] distributions particle distribution functions at a given point
	 * @return macroscopic density (sum of all distributions)
	 */
	double calculateDensity(const vector<double>& distributions) const {

		// calculate macroscopic density (rho)
		double rho = 0.0;
		for (size_t i = 0; i < m_q; i++){
			rho += distributions.at(i);
		}
		return rho;

	}

	/**
	 * @short calculate macroscopic velocity
	 * @param[in] distributions particle distribution functions at a given point
	 * @return macroscopic velocity
	 */
	numeric_vector calculateVelocity(const vector<double>& distributions) const {

		numeric_vector u(m_d);
		for (size_t i = 0; i < m_q; i++){
			// TODO efficient calculation of scalar*directions?
			Math::add_vector(u, Math::scalar_vector(distributions.at(i), m_directions.at(i)));
		}
		Math::scale_vector(1./calculateDensity(distributions), u);
		return u;
	}

	/**
	 * @short calculate macroscopic velocity; saves the double calculation of the density
	 * @note more efficient
	 * @param[in] distributions particle distribution functions at a given point
	 * @param[in] rho macroscopic density
	 * @param[out] u macroscopic velocity
	 */
	void calculateVelocity(const vector<double>& distributions, const double rho, numeric_vector& u) const {

		// assert
		assert(u.size() == m_d);
		assert(u(0) == 0.0);
		assert(u(m_d-1) == 0.0);

		for (size_t i = 0; i < m_q; i++){
			// TODO efficient calculation of scalar*directions?
			Math::add_vector(u, Math::scalar_vector(distributions.at(i), m_directions.at(i)));
		}
		Math::scale_vector(1./rho, u);
	}


//////////////////////////////
// EQUILIBRIUM DISTRIBUTION //
//////////////////////////////

	/** @short virtual function for the calculation of the equilibrium distribution
	 *  @param i index of the direction
	 *  @param u macroscopic velocity
	 *  @param rho macroscopic density
	 *  @return value of the equilibrium distribution
	 *  @note The calculation can surely be done more efficiently by passing different arguments,
	 *        e.g. u*u or u/(c^2)
	 */
	virtual double getEquilibriumDistribution(size_t i, const numeric_vector& u,
			const double rho = 1) const = 0;


	/** @short function for the calculation of all equilibrium distributions
	 *  @param[out] feq vector of all equality distributions, must have size Q
	 *  @param[in] u macroscopic velocity
	 *  @param[in] rho macroscopic density
	 *  @note The calculation can surely be done more efficiently by passing different arguments,
	 *        e.g. u*u or u/(c^2)
	 */
	virtual void getEquilibriumDistributions(vector<double>& feq, const numeric_vector& u,
			const double rho = 1) const;

	virtual double getSpeedOfSound() const = 0;
	virtual double getSpeedOfSoundSquare() const = 0;

};

} /* namespace natrium */
#endif /* BOLTZMANNMODEL_H_ */
