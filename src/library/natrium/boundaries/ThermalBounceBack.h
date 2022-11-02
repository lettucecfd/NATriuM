/*
 * LinearBoundaryRhoU.h
 *
 *  Created on: 29.10.2022
 *      Author: Dominik Wilde
 */

#ifndef LIBRARY_NATRIUM_PROBLEMDESCRIPTION_THERMALBOUNCEBACK_H_
#define LIBRARY_NATRIUM_PROBLEMDESCRIPTION_THERMALBOUNCEBACK_H_

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
class ThermalBounceBack: public Boundary<dim> {
public:
	/** @short This constructor assigns the Boundary condition with arbitrary density and velocity
	 *         to the boundary with the given boundary indicator.
	 *  @param[in] boundaryIndicator the boundary indicator that is assigned to the target boundary.
	 *  @param[in] boundaryVelocity A dealii::Function<dim> that defines the prescribed velocity at the boundary.
	 */
	ThermalBounceBack(size_t boundaryIndicator,
			boost::shared_ptr<dealii::Function<dim> > boundaryVelocity, double wallTemperature);

	/** @short This constructor assigns the Boundary condition with a constant fixed velocity and \f[ \rho = 1 \f]
	 *  to the boundary with the given boundary indicator.
	 *  @param[in] boundaryIndicator the boundary indicator that is assigned to the target boundary.
	 *  @param[in] velocity Constant velocity vector at the boundary.
	 */
	ThermalBounceBack(size_t boundaryIndicator,
			const dealii::Vector<double>& velocity, double wallTemperature);

	ThermalBounceBack(size_t boundaryIndicator,
			const dealii::Tensor<1,dim>& velocity, double wallTemperature);

	/// destructor
	virtual ~ThermalBounceBack();


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
		return false;
	}

	/** @short is the boundary set up to work with semi-Lagrangian streaming
	 */
	virtual bool isSLSupported() const{
		return true;
	}



	/////////////////////////////////////////////////
	// SL ///////////////////////////////////////////
	/////////////////////////////////////////////////

	virtual void calculateBoundaryValues(
			FEBoundaryValues<dim>& fe_boundary_values, size_t q_point,
			const LagrangianPathDestination& destination, double eps,
			double t);


	virtual BoundaryFlags getUpdateFlags() const {
		BoundaryFlags flags = boundary_rho;
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

private:
    double m_wallTemperature;

};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_PROBLEMDESCRIPTION_THERMALBOUNCEBACK_H_ */
