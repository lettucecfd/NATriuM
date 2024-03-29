/*
 * LinearBoundaryRhoU.h
 *
 *  Created on: 08.12.2015
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_PROBLEMDESCRIPTION_VELOCITYNEQBOUNCEBACK_H_
#define LIBRARY_NATRIUM_PROBLEMDESCRIPTION_VELOCITYNEQBOUNCEBACK_H_

#include "BoundaryTools.h"
#include "Boundary.h"

namespace natrium {

/*
 *			@short Boundary with fixed velocity (u) and density (\f[ \rho \f]).
 *			This boundary condition bounces back the nonequilibrium parts of the distribution
 *			function, where the equilibrium distribution is constructed with the target values.
 * 			For outgoing particle distribution functions the fluxes are set to 0.
 * 		  	For incoming particle distributions fluxes are set to
 * 		  	\f[ f_{\alpha} - f^{+}_{\alpha} = f_{\alpha} - f_{\alpha^{*}} - 2w_{\alpha} \rho_{0} (e_{\alpha}\cdot u_{b})/c^{2}_{s}\f]
 * 		  	This formula shows that we have to couple only opposite distribution functions (COUPLE_ONLY_OPPOSITE_DISTRIBUTIONS)
 * 		  	at single points (COUPLE_ONLY_SINGLE_POINTS - as there are no gradients to be computed)
 * 		  	@note: The boundary has been proposed by Ladd in the standard LBM and has also been used by Min and Lee (2011)
 * 		  	when proposing the SEDG-LBM. It has been shown practically that it  supports the exponential convergence of the scheme.
 */
template<size_t dim>
class VelocityNeqBounceBack: public Boundary<dim> {
public:
	/** @short This constructor assigns the Boundary condition with arbitrary density and velocity
	 *         to the boundary with the given boundary indicator.
	 *  @param[in] boundaryIndicator the boundary indicator that is assigned to the target boundary.
	 *  @param[in] boundaryVelocity A dealii::Function<dim> that defines the prescribed velocity at the boundary.
	 */
	VelocityNeqBounceBack(size_t boundaryIndicator,
			boost::shared_ptr<dealii::Function<dim> > boundaryVelocity);

	/** @short This constructor assigns the Boundary condition with a constant fixed velocity and \f[ \rho = 1 \f]
	 *  to the boundary with the given boundary indicator.
	 *  @param[in] boundaryIndicator the boundary indicator that is assigned to the target boundary.
	 *  @param[in] velocity Constant velocity vector at the boundary.
	 */
	VelocityNeqBounceBack(size_t boundaryIndicator,
			const dealii::Vector<double>& velocity);

	VelocityNeqBounceBack(size_t boundaryIndicator,
			const dealii::Tensor<1,dim>& velocity);

	/// destructor
	virtual ~VelocityNeqBounceBack();


	/////////////////////////////////////////////////
	// FLAGS ////////////////////////////////////////
	/////////////////////////////////////////////////

	/** @short is the boundary a periodic boundary ?
	 */
	virtual bool isPeriodic() const{
		return false;
	}

	/** @short is the boundary a linear flux boundary as in SEDG-LBM (affine linear in the distributions f)?
	 */
	virtual bool isDGSupported() const{
		return true;
	}

	/** @short is the boundary set up to work with semi-Lagrangian streaming
	 */
	virtual bool isSLSupported() const{
		return true;
	}


	/////////////////////////////////////////////////
	// SEDG /////////////////////////////////////////
	/////////////////////////////////////////////////
	/**
	 * @short This function defines the actual boundary condition. It calculates and assembles the fluxes at the
	 * boundary. Note that this function has to be called only once at the beginning of a simulation
	 * (or at reassembly, e.g. when the grid has changed).
	 * While the simulation runs, there has to be no separate treatment of the boundary, as it is a part
	 * of the linear ODE, already.
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
			bool useCentralFlux = false) ;

	virtual BoundaryTools::DistributionCouplingAtBoundary getDistributionCoupling() const {
		return  BoundaryTools::COUPLE_ONLY_OPPOSITE_DISTRIBUTIONS;
		}

	virtual BoundaryTools::PointCouplingAtBoundary getPointCoupling() const {
		return BoundaryTools::COUPLE_ONLY_SINGLE_POINTS;

	}

	void addToSparsityPattern(
			dealii::TrilinosWrappers::SparsityPattern& cSparse,
			const dealii::DoFHandler<dim>& doFHandler) const {
		BoundaryTools::CoupleDoFsAtBoundary<dim>(cSparse,
				doFHandler, Boundary<dim>::getBoundaryIndicator(), getPointCoupling());
	}


	/////////////////////////////////////////////////
	// SL ///////////////////////////////////////////
	/////////////////////////////////////////////////

	virtual void calculateBoundaryValues(
			FEBoundaryValues<dim>& fe_boundary_values, size_t q_point,
			const LagrangianPathDestination& destination, double eps,
			double t);


	virtual BoundaryFlags getUpdateFlags() const {
		BoundaryFlags flags = only_distributions;
		return flags;
	}


	/**
	 * @short Calculates outgoing distribution from incoming distributions.
	 * @param[in/out] boundary_hit The boundary hit instance that contains all information about the boundary hit.
	 * @param[in] stencil the stencil (e.g. a D2Q9 instance)
	 * @param[in] time_of_next_step the physical time at the next time step (is required here to define time-dependent boundary conditions)
	 * @note This function is used by the semi-Lagrangian advection solver. Before it is called on a
	 *       BoundaryHit instance, the BoundaryHit instance must have the right incoming directions (usually filled
	 *       by makeIncomingDirections()) and the right references in fIn (has to be filled by hand -- by the
	 *       semi-Lagrangian advection solver).
	 */
/*	virtual void calculate(BoundaryHit<dim>& boundary_hit,
			const Stencil& stencil, double time_of_next_step,
			SemiLagrangianVectorAccess& f) const {

		const numeric_vector ea = stencil.getDirection(
				boundary_hit.out.getAlpha());
		numeric_vector velocity(dim);
		LinearBoundary<dim>::getBoundaryVelocity()->set_time(
				time_of_next_step - boundary_hit.time_shift);
		LinearBoundary<dim>::getBoundaryVelocity()->vector_value(
				boundary_hit.coordinates, velocity);
		f[boundary_hit.out] = f(boundary_hit.in.at(0))
				+ 2.0 * stencil.getWeight(boundary_hit.out.getAlpha()) * 1.0
						* (ea * velocity) / stencil.getSpeedOfSoundSquare();
	}
*/
	/**
	 * @short Resizes boundary_hit.incomingDirections and fills it in.
	 * @param[in/out] boundary_hit The boundary hit instance that contains all information about the boundary hit.
	 * @param[in] stencil the stencil (e.g. a D2Q9 instance)
	 * @note This function is used by the semi-Lagrangian advection solver
	 */
/*	virtual void makeIncomingDirections(BoundaryHit<dim>& boundary_hit,
			const Stencil& stencil) const {
		assert (boundary_hit.incomingDirections.empty());
		boundary_hit.incomingDirections.push_back(stencil.getIndexOfOppositeDirection(
						boundary_hit.out.getAlpha()));
	}
	*/

};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_PROBLEMDESCRIPTION_VELOCITYNEQBOUNCEBACK_H_ */
