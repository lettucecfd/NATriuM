/*
 * DirichletBoundaryRhoU.h
 *
 *  Created on: 08.12.2015
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_PROBLEMDESCRIPTION_DIRICHLETBOUNDARYRHOU_H_
#define LIBRARY_NATRIUM_PROBLEMDESCRIPTION_DIRICHLETBOUNDARYRHOU_H_

#include "DirichletBoundary.h"

namespace natrium {

/*
 *			@short Boundary with fixed velocity (u) and density (\f[ \rho \f]).
 *			This boundary condition bounces back the nonequilibrium parts of the distribution
 *			function, where the equilibrium distribution is constructed with the target values.
 * 			For outgoing particle distribution functions the fluxes are set to 0.
 * 		  	For incoming particle distributions fluxes are set to
 * 		  	\f[ f_{\alpha} - f^{+}_{\alpha} = f_{\alpha} - f_{\alpha^{*}} - 2w_{\alpha} \rho_{0} (e_{\alpha}\cdot u_{b})/c^{2}_{s}\f]
 * 		  	@note: The boundary has been proposed by Ladd in the standard LBM and has also been used by Min and Lee (2011)
 * 		  	when proposing the SEDG-LBM. It has been shown practically that it  supports the exponential convergence of the scheme.
 */
template<size_t dim>
class DirichletBoundaryRhoU: public DirichletBoundary<dim> {
public:
	/** @short This constructor assigns the Boundary condition with arbitrary density and velocity
	 *         to the boundary with the given boundary indicator.
	 *  @param[in] boundaryIndicator the boundary indicator that is assigned to the target boundary.
	 *  @param[in] boundaryDensity A dealii::Function<dim> that defines the prescribed density at the boundary.
	 *  @param[in] boundaryVelocity A dealii::Function<dim> that defines the prescribed velocity at the boundary.
	 */
	DirichletBoundaryRhoU(size_t boundaryIndicator,
			boost::shared_ptr<dealii::Function<dim> > boundaryDensity,
			boost::shared_ptr<dealii::Function<dim> > boundaryVelocity);

	/** @short This constructor assigns the Boundary condition with a constant fixed velocity and \f[ \rho = 1 \f]
	 *  to the boundary with the given boundary indicator.
	 *  @param[in] boundaryIndicator the boundary indicator that is assigned to the target boundary.
	 *  @param[in] velocity Constant velocity vector at the boundary.
	 */
	DirichletBoundaryRhoU(size_t boundaryIndicator,
			const dealii::Vector<double>& velocity);

	/// destructor
	virtual ~DirichletBoundaryRhoU();

	/**
	 * @short In order to include the boundary conditions into the integration
	 * we have to expand the sparsity pattern. In concrete, we have to couple blocks at the boundary
	 * that belong to different distribution functions \f[ f_{\alpha} \f] and \f[ f_{\beta} \f].
	 * Here, we have to couple only the blocks of opposite distributions, as both density and velocity
	 * are prescribed.
	 */
	virtual void addToSparsityPattern(
			dealii::TrilinosWrappers::SparsityPattern& cSparse,
			const dealii::DoFHandler<dim>& doFHandler, const Stencil&) const;

	/**
	 * @short This function defines the actual boundary condition. It calculates and assembles the fluxes at the
	 * boundary. Note that this function has to be called only once at the beginning of a simulation
	 * (or at reassembly, e.g. when the grid has changed).
	 * While the simulation runs, there has to be no separate treatment of the boundary, as it is a part
	 * of the linear ODE, already.
	 */
	virtual void assembleBoundary(size_t alpha,
			const typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
			size_t faceNumber, dealii::FEFaceValues<dim>& feFaceValues,
			const Stencil& stencil,
			const std::map<size_t, size_t>& q_index_to_facedof,
			const vector<double> & inverseLocalMassMatrix,
			distributed_sparse_block_matrix& systemMatrix,
			distributed_block_vector& systemVector,
			bool useCentralFlux = false) const;

};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_PROBLEMDESCRIPTION_DIRICHLETBOUNDARYRHOU_H_ */
