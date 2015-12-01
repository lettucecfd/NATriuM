/**
 * @file MinLeeBoundary.h
 * @short Description of a boundary as described by Min and Lee
 * @date 26.03.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef MINLEEBOUNDARY_H_
#define MINLEEBOUNDARY_H_

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include "deal.II/lac/trilinos_sparsity_pattern.h"
#include "deal.II/dofs/dof_handler.h"
#include "deal.II/base/function.h"

#include "Boundary.h"
#include "../stencils/Stencil.h"
#include "../utilities/BasicNames.h"

namespace natrium {

template<size_t dim>
class BoundaryDensity: public dealii::Function<dim> {
public:
	BoundaryDensity() {
	}
	;
	virtual ~BoundaryDensity() {
	}
	;
	virtual double value(const dealii::Point<dim> &,
			const unsigned int  = 0) const {
		return 1;
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
 * @short 	The boundary described by Min and Lee.
 * 			For outgoing particle distribution functions the fluxes are set to 0
 * 		  	For incoming particle distributions fluxes are set to
 * 		  	\f[ f_{\alpha} - f^{+}_{\alpha} = f_{\alpha} - f_{\alpha^{*}} - 2w_{\alpha} \rho_{0} (e_{\alpha}\cdot u_{b})/c^{2}_{s}\f]
 *
 */
template<size_t dim> class MinLeeBoundary: public Boundary<dim> {
private:

	size_t m_boundaryIndicator;

	boost::shared_ptr<dealii::Function<dim> > m_boundaryDensity;

	boost::shared_ptr<dealii::Function<dim> > m_boundaryVelocity;

public:

	/// constructor
	MinLeeBoundary(size_t boundaryIndicator,
			boost::shared_ptr<dealii::Function<dim> > boundaryDensity,
			boost::shared_ptr<dealii::Function<dim> > boundaryVelocity);

	/// constructor
	MinLeeBoundary(size_t boundaryIndicator,
			const dealii::Vector<double>& velocity);

	/// destructor
	virtual ~MinLeeBoundary() {
	}
	;

	/**
	 * @short modify sparsity pattern so that the fluxes over periodic boundary can be incorporated
	 * @param cSparse the block-sparsity pattern
	 */

	void addToSparsityPattern(
#ifdef WITH_TRILINOS_MPI
			dealii::TrilinosWrappers::SparsityPattern& cSparse,
#else

			dealii::DynamicSparsityPattern& cSparse,
#endif
			const dealii::DoFHandler<dim>& doFHandler,
			const Stencil& ) const;

	// assemble Min-Lee-Type boundary
	void assembleBoundary(size_t alpha,
			const typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
			size_t faceNumber, dealii::FEFaceValues<dim>& feFaceValues,
			const Stencil& stencil,
			const std::map<size_t, size_t>& q_index_to_facedof,
			const vector<double> & inverseLocalMassMatrix,
			distributed_sparse_block_matrix& systemMatrix,
			distributed_block_vector& systemVector,
			bool useCentralFlux = false) const;

	size_t getBoundaryIndicator() const {
		return m_boundaryIndicator;
	}
};

} /* namespace natrium */

#endif /* MINLEEBOUNDARY_H_ */
