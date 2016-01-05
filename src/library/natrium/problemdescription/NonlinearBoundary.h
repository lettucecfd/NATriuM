/**
 * @file NonlinearBoundary.h
 * @short Description of a boundary that cannot be described in the linear framework of Min and Lee
 * @date 15.12.2015
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef NONLINEAR_BOUNDARY_H_
#define NONLINEAR_BOUNDARY_H_

#include "deal.II/lac/dynamic_sparsity_pattern.h"
#include "deal.II/lac/trilinos_sparsity_pattern.h"
#include "deal.II/base/function.h"
#include "deal.II/fe/fe_update_flags.h"

#include "Boundary.h"
#include "BoundaryTools.h"
#include "../advection/AdvectionOperator.h"
#include "../solver/DistributionFunctions.h"
#include "../stencils/Stencil.h"
#include "../utilities/BasicNames.h"

namespace natrium {


template<size_t dim>
class BoundaryDensity: public dealii::Function<dim> {
private:
	double m_density;
public:
	BoundaryDensity(double rho = 1) {
		m_density = rho;
	}
	;
	virtual ~BoundaryDensity() {
	}
	;
	virtual double value(const dealii::Point<dim> &,
			const unsigned int  = 0) const {
		return m_density;
	}
};
template<size_t dim>
class BoundaryVelocity: public dealii::Function<dim> {
private:
	dealii::Vector<double> m_Velocity;
public:
	BoundaryVelocity(const dealii::Vector<double>& velocity) :
			m_Velocity(velocity) {
	}
	virtual ~BoundaryVelocity() {
	}
	;
	virtual void vector_value(const dealii::Point<dim> &,
			dealii::Vector<double> &values) const {
		values = m_Velocity;
	}
};

/**
 * @short 	Abstract class to describe Linear boundary conditions.
 * 			The virtual function to be overriden is assembleBoundary. Moreover, the DoF couplings at the
 * 			boundary have to be defined (see Documentation of the constructors).
 */
template<size_t dim> class NonlinearBoundary: public Boundary<dim> {
private:

	size_t m_boundaryIndicator;
	dealii::UpdateFlags m_updateFlags;
	boost::shared_ptr<AdvectionOperator<dim> > m_advectionOperator;
	boost::shared_ptr<Stencil> m_stencil;
	distributed_vector const * m_rho;
	vector<distributed_vector> const * m_u;
	DistributionFunctions const * m_f;
	distributed_block_vector* m_boundaryVector;

public:

	/**
	 * Constructor
	 *  @param[in] boundaryIndicator the boundary indicator that is assigned to the target boundary.
	 *  @param[in] boundaryDensity A dealii::Function<dim> that defines the prescribed density at the boundary.
	 *  @param[in] boundaryVelocity A dealii::Function<dim> that defines the prescribed velocity at the boundary.
	 */
	NonlinearBoundary(size_t boundaryIndicator, const dealii::UpdateFlags update_flags = dealii::update_values
			| dealii::update_JxW_values | dealii::update_normal_vectors);


	/// destructor
	virtual ~NonlinearBoundary() {
	}
	;

	virtual void updateNonlinearBoundaryValues() const = 0;

	void initialize(boost::shared_ptr<AdvectionOperator<dim> > advection_operator, boost::shared_ptr<Stencil> stencil, const distributed_vector* rho,
			const vector<distributed_vector>* u, const DistributionFunctions* f,
			distributed_block_vector* boundary_vector);

	size_t getBoundaryIndicator() const {
		return m_boundaryIndicator;
	}

	/** @short is the boundary a dirichlet boundary ?
	 */
	virtual bool isLinear() const {
		return false;
	}

	const AdvectionOperator<dim>& getAdvectionOperator() const {
		return *m_advectionOperator;
	}

	distributed_block_vector& getBoundaryVector() const {
		return *m_boundaryVector;
	}

	const DistributionFunctions& getF() const {
		return *m_f;
	}

	const distributed_vector& getRho() const {
		return *m_rho;
	}

	const vector<distributed_vector>& getU() const {
		return *m_u;
	}

	const dealii::UpdateFlags& getUpdateFlags() const {
		return m_updateFlags;
	}

	const Stencil& getStencil() const {
		return *m_stencil;
	}

	void resetBoundaryVector() {
		*m_boundaryVector = m_advectionOperator->getSystemVector();
	}
};

} /* namespace natrium */

#endif /* LinearBoundary_H_ */
