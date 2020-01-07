/**
 * @file Stencil.h
 * @short Abstract class for the description of a boltzmann model
 * @date 04.06.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef STENCIL_H_
#define STENCIL_H_

#include "assert.h"

#include "../utilities/Math.h"
#include "../utilities/BasicNames.h"

namespace natrium {

/**
 * @short Enum type for the difference stencil
 */
enum StencilType {

	/// D2Q9 stencil
	Stencil_D2Q9,
	// D3Q19 stencil
	Stencil_D3Q19,
	// D3Q15 stencil
	Stencil_D3Q15,
	// D3Q27 stencil
	Stencil_D3Q27,
	/// D2Q25 stencil
	Stencil_D2Q25,
	Stencil_D2Q25H,
	Stencil_D3Q27,
	// D3Q13 stencil
	Stencil_D3Q13,
	// D3Q21 non-cubic stencil
	Stencil_D3Q21,
	// Crystallographic RD3Q27 stencil
	Stencil_RD3Q27

};

/**
 * @short Abstract class for the description of a boltzmann model, e.g. D2Q9
 */
class Stencil {

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

	// transformation matrix to project the distribution functions
	// onto an appropriate Q-dimensional subspace of the moment space
	numeric_matrix m_momentBasis;
	numeric_matrix m_inverseMomentBasis;


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
	Stencil(size_t d, size_t q, const vector<numeric_vector>& directions,
			const vector<double>& weights, StencilType stencilType, const numeric_matrix& moment_basis);

	/// destructor
	virtual ~Stencil();

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
	const StencilType& getStencilType() const {
		return m_stencilType;
	}

	virtual size_t getIndexOfOppositeDirection(size_t index) const = 0;


	virtual double getSpeedOfSound() const = 0;
	virtual double getSpeedOfSoundSquare() const = 0;
	virtual double getMaxParticleVelocityMagnitude() const = 0;

	void getMomentBasis(numeric_matrix& m) const{
		m = m_momentBasis;
	}
	void getInverseMomentBasis(numeric_matrix& m_inv) const {
		m_inv = m_inverseMomentBasis;
	}

	virtual double getScaling() const {
		return 1.0;
	}

};

} /* namespace natrium */
#endif /* STENCIL_H_ */
