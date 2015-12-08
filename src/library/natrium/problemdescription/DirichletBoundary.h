/**
 * @file DirichletBoundary.h
 * @short Description of a boundary as described by Min and Lee
 * @date 26.03.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef DIRICHLETBOUNDARY_H_
#define DIRICHLETBOUNDARY_H_

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include "deal.II/lac/trilinos_sparsity_pattern.h"
#include "deal.II/dofs/dof_handler.h"
#include "deal.II/base/function.h"

#include "Boundary.h"
#include "BoundaryTools.h"
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
 * @short 	Abstract class to describe Dirichlet boundary conditions.
 * 			The virtual function to be overriden is assembleBoundary. Moreover, the DoF couplings at the
 * 			boundary have to be defined (see Documentation of the constructors).
 */
template<size_t dim> class DirichletBoundary: public Boundary<dim> {
private:

	size_t m_boundaryIndicator;

	boost::shared_ptr<dealii::Function<dim> > m_boundaryDensity;

	boost::shared_ptr<dealii::Function<dim> > m_boundaryVelocity;

public:

	const BoundaryTools::DistributionCouplingAtBoundary m_distributionCoupling;

	const BoundaryTools::PointCouplingAtBoundary m_pointCoupling;

	/**
	 * Constructor
	 *  @param[in] boundaryIndicator the boundary indicator that is assigned to the target boundary.
	 *  @param[in] boundaryDensity A dealii::Function<dim> that defines the prescribed density at the boundary.
	 *  @param[in] boundaryVelocity A dealii::Function<dim> that defines the prescribed velocity at the boundary.
	 *  @param[in] distribution_coupling Indicates which distributions are coupled at the boundary. There are two possibilities:
	 *  		   -# \[ f_{\alpha} \f] depends only on the opposite distribution function  \[ f_{\alpha^{\ast}} \f]. Then this
	 *  		   parameter has to be COUPLE_ONLY_OPPOSITE_DISTRIBUTIONS.
	 *  		   -#  \[ f_{\alpha} \f] depens on all distribution functions (e.g. when the density or velocity is
	 *  		   used for the boundary condition (the calculated on, not the prescribed one). Then, this parameter has to be
	 *  		   COUPLE_ALL_DISTRIBUTIONS.
	 * @param[in] point_coupling Indicates which DoFs are coupled at the boundary. There are two cases:
	 * 			  -# \[ f_{\alpha} \f] depends only on the distribution functions at a given point (and possibly some other data, but no
	 * 			  distribution functions at other points). Then, this parameter has to be COUPLE_ONLY_SINGLE_POINTS.
	 * 			  -#  \[ f_{\alpha} \f] depends on the distribution functions at other points at the face (e.g. when gradients
	 * 			  are computed.) Then, this parameter has to be COUPLE_WHOLE_FACE.
	 */
	DirichletBoundary(size_t boundaryIndicator,
			boost::shared_ptr<dealii::Function<dim> > boundaryDensity,
			boost::shared_ptr<dealii::Function<dim> > boundaryVelocity,
			BoundaryTools::DistributionCouplingAtBoundary distribution_coupling,
			BoundaryTools::PointCouplingAtBoundary point_coupling);

	/** @short This constructor assigns the Boundary condition with a constant fixed velocity and \f[ \rho = 1 \f]
	 *  to the boundary with the given boundary indicator.
	 *  @param[in] boundaryIndicator the boundary indicator that is assigned to the target boundary.
	 *  @param[in] velocity Constant velocity vector at the boundary.
	 *  @param[in] distribution_coupling See previous constructor.
	 *  @param[in] point_coupling See previous constructor.
	 */
	DirichletBoundary(size_t boundaryIndicator,
			const dealii::Vector<double>& velocity,
			BoundaryTools::DistributionCouplingAtBoundary distribution_coupling,
			BoundaryTools::PointCouplingAtBoundary point_coupling);

	/// destructor
	virtual ~DirichletBoundary() {
	}
	;

	/**
	 * @short In order to include the boundary conditions into the integration
	 * we have to expand the sparsity pattern. In concrete, we have to couple blocks at the boundary
	 * that belong to different distribution functions \f[ f_{\alpha} \f] and \f[ f_{\beta} \f].
	 * The concrete coupling depends on m_pointCoupling. If this is
	 * - COUPLE_ONLY_SINGLE_POINTS: only the DoFs belonging to the same point are coupled
	 * - COUPLE_WHOLE_FACE: all DoFs at the face are coupled, e.g. if \f[ \partial_{x}f_{\beta} \f]
	 *   is required to calculate \f[ f_{\alpha} \f] at the boundary.
	 * @param[in/out] doFHandler the DoFHandler
	 *
	 */
	virtual void addToSparsityPattern(
			dealii::TrilinosWrappers::SparsityPattern& cSparse,
			const dealii::DoFHandler<dim>& doFHandler) const;

	/**
	 * @short Pure virtual function for the concrete assembly of the boundary.
	 *        The child classes have to override this function for each specific boundary.
	 * @param[in] alpha The index of the distribution function \f f_{\alpha} \f] for which the boundary is constructed.
	 * @param[in] cell The current cell (assembly is done cell-by-cell)
	 * @param[in] faceNumber The current face number.
	 * @param[in] feFaceValues The dealii::FEFaceValues<dim> object, which has to be reinitialized right at the beginning
	 *            of the assembleBoundary function (it can be empty, but is passed here as a function parameter to avoid
	 *            construction of a new feFaceValues object).
	 * @param[in] stencil The DQ stencil, which is required here to get the particle velocities.
	 * @param[in] q_index_to_facedof Maps quadrature nodes to DoFs (only possible for GLL nodes)
	 * @param[in] inverseLocalMassMatrix \f[ M^{-1} \f] is multiplied with the obtained boundary matrix at the end of the assembly.
	 * @param[in/out] systemMatrix The global system matrix which is assembled to.
	 * @param[in/out] systemVector The global system vector which is assembled to.
	 * @param[in] useCentralFlux indicates whether to use a central instead of a Lax-Friedrichs flux. Should not be used,
	 *            has not been tested thoroughly and yields bad results, usually.
	 */
	virtual void assembleBoundary(size_t alpha,
			const typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
			size_t faceNumber, dealii::FEFaceValues<dim>& feFaceValues,
			const Stencil& stencil,
			const std::map<size_t, size_t>& q_index_to_facedof,
			const vector<double> & inverseLocalMassMatrix,
			distributed_sparse_block_matrix& systemMatrix,
			distributed_block_vector& systemVector,
			bool useCentralFlux = false) const = 0;

	size_t getBoundaryIndicator() const {
		return m_boundaryIndicator;
	}

	const boost::shared_ptr<dealii::Function<dim> >& getBoundaryDensity() const {
		return m_boundaryDensity;
	}

	const boost::shared_ptr<dealii::Function<dim> >& getBoundaryVelocity() const {
		return m_boundaryVelocity;
	}
};

} /* namespace natrium */

#endif /* DirichletBoundary_H_ */
