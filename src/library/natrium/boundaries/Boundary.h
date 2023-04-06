/**
 * @file Boundary.h
 * @short Abstract class for Description of a Boundary object
 * @date 24.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef NATRIUM_BOUNDARIES_BOUNDARY_H_
#define NATRIUM_BOUNDARIES_BOUNDARY_H_

#include "deal.II/grid/tria.h"
#include "deal.II/dofs/dof_handler.h"
#include "deal.II/fe/fe_values.h"
#include "deal.II/lac/dynamic_sparsity_pattern.h"
#include "deal.II/lac/trilinos_sparsity_pattern.h"
#include "deal.II/dofs/dof_handler.h"
#include "deal.II/base/function.h"

#include "Boundary.h"
#include "BoundaryFlags.h"
#include "FEBoundaryValues.h"
#include "BoundaryTools.h"
#include "../stencils/Stencil.h"
#include "../utilities/BasicNames.h"

namespace natrium {


enum BoundaryName {
	BOUNDARY_NOT_SET,
	PERIODIC_BOUNDARY,
	// the next three are combined in a single class VelocityNeqBounceBack
	// however, zero and constant velocity BB allow for more efficient implementations
	ZERO_VELOCITY_NEQ_BB,
	CONSTANT_VELOCITY_NEQ_BB, // constant means constant in time
	NONCONSTANT_VELOCITY_NEQ_BB,
	//
	FIRST_ORDER_BB,
	PRESSURE_EQUILIBRIUM_BOUNDARY,
	VELOCITY_EQUILIBRIUM_BOUNDARY,
	DO_NOTHING_BC,
    THERMAL_BB
};

inline bool is_velocity_neq_bb(BoundaryName bn){
	return ( (bn == ZERO_VELOCITY_NEQ_BB)
			or (bn == CONSTANT_VELOCITY_NEQ_BB)
			or (bn == NONCONSTANT_VELOCITY_NEQ_BB)
            or (bn == THERMAL_BB) );
}

    inline bool is_do_nothing_bb(BoundaryName bn){
        return bn == DO_NOTHING_BC;
    }

/**
 * @short  Abstract class for the description of boundaries.
 *         Base class for all boundaries.
 *         NATriuM supports two different advection schemes, SemiLagrangian
 *         and SEDG.
 *         In order for a boundary to support the SemiLagrangian streaming,
 *         it has to override the functions getUpdateFlags() and calculateBoundaryValues(...).
 *         In order for a boundary to support the SEDG streaming,
 *         it has to override the functions getDistributionCoupling(), getPointCopuling(),
 *         addToSparsityPattern(...), and assembleBoundary(...).
 *     	   Periodic boundaries present a special case. They do not have to
 *     	   override either of these functions.
 * @tparam dim The dimension of the boundary is the dimension of the domain -1 (
 * 	       e.g. 2-dim meshes have 1-dim boundary)
 */
template<size_t dim> class Boundary {
private:
	size_t m_boundaryIndicator;
	BoundaryName m_boundaryName;
	PrescribedBoundaryValues<dim> m_boundaryValues;

protected:
	/** @short subclasses can change their boundary name
	 */
	void setBoundaryName(BoundaryName boundary_name){
		m_boundaryName = boundary_name;
	}

public:

	/// constructor
	Boundary(size_t boundary_indicator, BoundaryName boundary_name,
			const PrescribedBoundaryValues<dim>& values):
		m_boundaryIndicator(boundary_indicator),
		m_boundaryName(boundary_name),
		m_boundaryValues(values){
	};

	/// destructor
	virtual ~Boundary(){};

	/** @short get the boundary id
	 */
	size_t getBoundaryIndicator()  const {
		return m_boundaryIndicator;
	}

	/** @short get the boundary name
	 */
	BoundaryName getBoundaryName() const {
		return m_boundaryName;
	}

	PrescribedBoundaryValues<dim>& getBoundaryValues()  {
		return m_boundaryValues;
	}

    BoundaryFlags& getPrescribedQuantities() const {
		return m_boundaryValues.getPrescribedValues();
	}

	/////////////////////////////////////////////////
	// FLAGS ////////////////////////////////////////
	/////////////////////////////////////////////////

	/** @short is the boundary a periodic boundary ?
	 */
	virtual bool isPeriodic() const = 0;

	/** @short is the boundary a linear flux boundary as in SEDG-LBM (affine linear in the distributions f)?
	 */
	virtual bool isDGSupported() const = 0;

	/** @short is the boundary set up to work with semi-Lagrangian streaming
	 */
	virtual bool isSLSupported() const = 0;



	/////////////////////////////////////////////////
	// SEDG /////////////////////////////////////////
	/////////////////////////////////////////////////


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
			const dealii::DoFHandler<dim>& doFHandler) const {
		throw NotImplementedException("To work with SEDG, this boundary must override addToSparsityPattern(...).");
	}

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
			bool useCentralFlux = false) {
		std::stringstream err;
		err << " The assembleBoundary() function is not implemented for this type of boundary "
				"but is required if you want to use the boundary with an SEDG solver (boundary id:"
				<< endl;
		throw NotImplementedException(err.str());
	}

	virtual BoundaryTools::DistributionCouplingAtBoundary getDistributionCoupling() const {
		throw NotImplementedException("To work with SEDG, this boundary must override getDistributionCoupling().");
	}

	virtual BoundaryTools::PointCouplingAtBoundary getPointCoupling() const {
		throw NotImplementedException("To work with SEDG, this boundary must override getPointCoupling().");

	}

	/////////////////////////////////////////////////
	// SEMI-LAGRANGIAN //////////////////////////////
	/////////////////////////////////////////////////

	/**
	 * @short Get update flags. This pure virtual function has to be overriden by derived classes.
	 * @return the flags specifying the required values that have to be calculated from the flow field at t-dt,
	 * 			e.g. only_distributions, boundary_rho, boundary_drho_dt, ...
	 */
	virtual BoundaryFlags getUpdateFlags() const  {
		throw NATriuMException("SL streaming can only be used when getUpdateFlags is "
				"overriden in Boundary definition.");
	}


	/**
	 * @short calculate boundary values. This pure virtual function has to be overriden by derived classes.
	 * @param fe_boundary_values An instance of FEBoundaryValues, that stores the flow variables at t-dt (here only the distribution functions)
	 * @param q_point the local index of the boundary hit point, i.e. its position in the fe_boundary_values
	 * @param destination defines the degree of freedom and discrete direction that the calculated value is assigned to
	 * 			(usually defines a point close to the boundary)
	 * @param eps  Defines the point in time at which the boundary is hit.
	 * @param t global time. Only relevant for time-dependent boundary conditions.
	 */
	virtual void calculateBoundaryValues(FEBoundaryValues<dim>& fe_boundary_values,
				size_t q_point, const LagrangianPathDestination& destination,
				double eps, double t) {
		throw NATriuMException("SL streaming can only be used when calculateBoundaryValues is "
						"overriden in Boundary definition.");
	}


};



} /* namespace natrium */

#endif /* NATRIUM_BOUNDARIES_BOUNDARY_H_ */
