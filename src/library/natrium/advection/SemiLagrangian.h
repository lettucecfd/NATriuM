/**
 * @file SemiLagrangian.h
 * @short Semi-Lagrangian advection operator
 * @date 29.04.16
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef SEMILAGRANGIAN_H_
#define SEMILAGRANGIAN_H_

#include <map>
#include <array>

#include "deal.II/grid/tria.h"
#include "deal.II/fe/fe_dgq.h"
#include "deal.II/dofs/dof_handler.h"
#include "deal.II/dofs/dof_tools.h"
#include "deal.II/lac/sparse_matrix.h"
#include "deal.II/fe/fe_update_flags.h"
#include "deal.II/fe/mapping_cartesian.h"
#include "deal.II/lac/block_sparsity_pattern.h"
#include "deal.II/base/quadrature_lib.h"
#include "deal.II/base/utilities.h"

#include "AdvectionOperator.h"
#include "AdvectionTools.h"
#include "SemiLagrangianTools.h"
#include "../boundaries/SemiLagrangianBoundaryHandler.h"
#include "../problemdescription/BoundaryCollection.h"
#include "../utilities/BasicNames.h"
#include "../utilities/NATriuMException.h"
#include "../utilities/Timing.h"
#include "../utilities/Logging.h"

namespace natrium {

/* forward declaration */
class Stencil;

/**
 * @short Exception class for SemiLagrangian advection operator
 */
class SemiLagrangianException: public NATriuMException {
private:
	std::string message;
public:
	SemiLagrangianException(const char *msg) :
			NATriuMException(msg), message(msg) {
	}
	SemiLagrangianException(const string& msg) :
			NATriuMException(msg), message(msg) {
	}
	~SemiLagrangianException() throw () {
	}
	const char *what() const throw () {
		return this->message.c_str();
	}
};

/** @short This class solves the linear advection equations by a semi-Lagrangian scheme
 * @tparam dim The dimension of the flow (2 or 3).
 */
template<size_t dim> class SemiLagrangian: public AdvectionOperator<dim> {

private:

	// quick access to data of base class (without this-> pointer)
	typedef AdvectionOperator<dim> Base;

	/// Sparsity Pattern of the sparse matrix
	/// The sparsity pattern is only used in the assembly
	/// storing only one line of blocks at a time significantly reduces the storage
	dealii::TrilinosWrappers::BlockSparsityPattern m_blockSparsityPattern;


	SemiLagrangianBoundaryHandler<dim> m_boundaryHandler;

	/**
	 * @short update the sparsity pattern of the system matrix // the sparse matrix
	 * @param the row index is used in assembly of the sparsity pattern (we want
	 *        to store only one line at a time to save memory)
	 */
	void fillSparseObject(bool sparsity_pattern = false, size_t row_index  = -1);

	/**
	 * @short update the sparsity pattern of the system matrix
	 */
	void updateSparsityPattern();

public:

	/// constructor
	/**
	 * @short Constructor
	 * @param[in] triangulation The global mesh.
	 * @param[in] orderOfFiniteElement The number of nodes element and dimension
	 * @param[in] stencil the DQ model
	 * @param[in] delta_t time step size; if delta_t = 0, the sparsity pattern is not updated during construction
	 */
	SemiLagrangian(ProblemDescription<dim>& problem,
			size_t orderOfFiniteElement, QuadratureName quad_name,
			SupportPointsName points_name, boost::shared_ptr<Stencil> stencil,
			double delta_t);

	SemiLagrangian(ProblemDescription<dim>& problem,
			size_t orderOfFiniteElement, boost::shared_ptr<Stencil> stencil,
			double delta_t);

	/*SemiLagrangian(boost::shared_ptr<ProblemDescription<dim> > problem, SupportPointsName quad_name,
	 size_t orderOfFiniteElement, boost::shared_ptr<Stencil> stencil,
	 double delta_t);*/

	/// destructor
	virtual ~SemiLagrangian() {
	}
	;

	/**
	 * @short Determines which face is crossed first, when moving from one point inside the cell to a point outside.
	 * @param[in] cell iterator to the active cell that contains the point p_inside
	 * @param[in] p_inside the point inside the cell
	 * @param[in] p_outside the point outside  the cell
	 * @param[out] p_boundary the point where the boundary is hit
	 * @param[out] lambda the parameter lambda that solves   p_boundary_unit = lambda * p_outside_unit + (1-lambda) * p_inside_unit
	 * @return face_id, if a face is crossed; -1, if no face is crossed (i.e. the second point is inside the cell)
	 * @note lambda is calculated for the unit cell. In general, p_boundary = lambda * p_outside + (1-lambda) * p_inside does not hold
	 * @note The current implementation does not do anything special at corner nodes. It prefers x over y over z faces. This may lead to problems later on.
	 */
	int faceCrossedFirst(
			const typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
			const dealii::Point<dim>& p_inside,
			const dealii::Point<dim>& p_outside, dealii::Point<dim>& p_boundary,
			double* lambda, size_t* child_id, bool without_mapping = false);

	/**
	 * @short This function does the same as faceCrossedFirst, but without using mapping functions
	 * @note The reason for this implementation is that the mapping function may not be invertible outside the cell
	 */
	int faceCrossedFirstWithoutMapping(
			const typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
			const dealii::Point<dim>& p_inside,
			const dealii::Point<dim>& p_outside, dealii::Point<dim>& p_boundary,
			double* lambda, size_t* child_id);

	/// function to (re-)assemble linear system
	virtual void reassemble();

	virtual void setupDoFs();

	virtual double stream(DistributionFunctions& f_old,
			DistributionFunctions& f, double t) {

		assert(&f_old != &f);
		f_old = f;
		Base::m_systemMatrix.vmult(f.getFStream(), f_old.getFStream());
		{
			TimerOutput::Scope timer_section(Timing::getTimer(), "Semi-Lagrangian Boundaries");
			m_boundaryHandler.apply(f, f_old, t);
		}
		return Base::m_deltaT;
	}

	virtual void applyBoundaryConditions(DistributionFunctions& f_old, DistributionFunctions& f, double t){
        m_boundaryHandler.apply(f, f_old, t);
    }

    virtual void applyBoundaryConditions(DistributionFunctions& f_old, DistributionFunctions& f, DistributionFunctions& g, double t){
        m_boundaryHandler.apply(f, f_old, g, t);
    }

	virtual void applyBoundaryConditionsToG(DistributionFunctions& f, DistributionFunctions& g, double t, const double gamma){
        m_boundaryHandler.applyToG(f, g, t, gamma);
    }

	virtual void setDeltaT(double deltaT) {
		Base::setDeltaT(deltaT);
		m_boundaryHandler.setTimeStep(deltaT);
		updateSparsityPattern();
	}

	/**
	 * @short get i-th neighbor of a cell, incorporating periodic boundaries
	 * @param cell An iterator pointing to cell
	 * @param i face index (=neighbor index)
	 * @return m_doFHandler->end(), if cell has no i-th neighbor (e.g. at solid boundary)
	 *
	 */
	typename dealii::DoFHandler<dim>::cell_iterator getNeighbor(
			const typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
			size_t i);

	/**
	 * @short fill the neighborhood list
	 * @param cell cell
	 * @param neighborhood the neighborhood object
	 * @note the neighborhood incorporates the current cell, all its neighbors,
	 * 		 and their respective neighbors; each cell has only pointer to it in the neighborhood.
	 */
	void getNeighborhood(
			typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
			Neighborhood<dim>& neighborhood, size_t n_shells = 1);

	/**
	 * @short recursively search a point in neighborhood, until is found
	 * @param p the point you search for
	 * @param cell The start cell of the recursive search
	 * @return A cell that contains the point p. If the point was not found, the cell pointer will point to DoFHandler.end()
	 */
	typename dealii::DoFHandler<dim>::active_cell_iterator recursivelySearchInNeighborhood(
			const dealii::Point<dim>& p,
			typename dealii::DoFHandler<dim>::active_cell_iterator& cell);

	virtual void setTimeIntegrator(
			boost::shared_ptr<
					TimeIntegrator<distributed_sparse_block_matrix,
							distributed_block_vector> >) {

	}

	const SemiLagrangianBoundaryHandler<dim>& getBoundaryHandler() const {
		return m_boundaryHandler;
	}

	/**
	 * @short Checks if a cell is already in the neighborhood list
	 * @param cell the cell
	 * @param neighborhood the neighborhood list
	 * @return true, if cell is already in the list
	 */
	bool isCellInNeighborhood(
			const typename dealii::DoFHandler<dim>::cell_accessor& cell,
			const Neighborhood<dim>& neighborhood) {
		for (size_t i = 0; i < neighborhood.size(); i++) {
			if (cell.id() == neighborhood.at(i)->id()) {
				return true;
			}
		}
		return false;
	}



	virtual const distributed_block_vector& getSystemVector() const {
		throw AdvectionSolverException("getSytemVector is not defined for SemiLagrangian streaming. Function to be removed."
				"as soon as AdvectionOperator::stream() works.");
	}

}
;

} /* namespace natrium */

#endif /* SEMILAGRANGIAN_H_ */
