/**
 * @file PeriodicBoundary.h
 * @short Description of a periodic boundary on a line.
 * @date 25.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef PERIODICBOUNDARY1D_H_
#define PERIODICBOUNDARY1D_H_

#include <exception>
#include <string>
#include <functional>
#include <map>

#include "deal.II/base/point.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"
#include "deal.II/dofs/dof_handler.h"

#include "Boundary.h"

namespace natrium {

/**
 * @short Exception for impossible periodic boundary,
 * e.g. when the interfaces don't have the same length.
 */
class PeriodicBoundaryNotPossible: public std::exception {
private:
	std::string m_message;
public:
	PeriodicBoundaryNotPossible() :
			m_message("Error in periodic boundary") {
	}
	PeriodicBoundaryNotPossible(const std::string& message) :
			m_message("Error in periodic boundary: " + message) {
	}
	virtual const char* what() const throw () {
		return m_message.c_str();
	}
	virtual ~PeriodicBoundaryNotPossible() throw () {
	}
};


/**
 * @short  A periodic boundary condition,
 * @note   First use in step-1 tutorial.
 */
template<size_t dim>
class PeriodicBoundary: public Boundary<dim> {
private:

	/// triangulation object
	shared_ptr<dealii::Triangulation<dim> > m_triangulation;

	/// Container for all cells that belong to this boundary
	/// stored as <accessor to cell, (accessor to opposite cell, boundary face at opposite cell) > /
	std::map<typename dealii::DoFHandler<dim>::active_cell_iterator, std::pair<typename dealii::DoFHandler<dim>::cell_iterator, size_t> > m_cells;

	/// boundary indicator of first interfacial line
	size_t m_boundaryIndicator1;

	/// boundary indicator of second interfacial line
	size_t m_boundaryIndicator2;


	//////////////////////////
	// ONLY RELEVANT FOR 2D //
	//////////////////////////
	/// start point of line 1
	dealii::Point<dim> m_beginLine1;

	/// end point of line 1
	dealii::Point<dim> m_endLine1;

	/// start point of line 2
	dealii::Point<dim> m_beginLine2;

	/// end point of line 2
	dealii::Point<dim> m_endLine2;


	/**
	 * @short Check if the two lines are OK (right positions, lengths, etc).
	 *        If they are anti-parallel, begin and end vector of the second are swapped.
	 *
	 * @throws PeriodicBoundaryNotPossible exception if not OK
	 */
	void checkInterfacePositions();

	/**
	 * @short Get the positions of the lines that define interfaces from the
	 *        defining Boundary indicators.
	 *
	 *  @param boundaryIndicator1 boundary indicator of interface line 1
	 *  @param boundaryIndicator2 boundary indicator of interface line 2
	 *  @param triangulation A (shared ptr to a) triangulation object (the mesh)
	 *  @param[out] beginLine1 start point of line 1
	 *  @param[out] endLine1 end point of line 1
	 *  @param[out] beginLine2 start point of line 2
	 *  @param[out] endLine2 end point of line 2
	 *
	 *  @throws PeriodicBoundaryNotPossible exception if not all points with same boundary indicator are on a line
	 */
	void getInterfacePositionsByBoundaryIndicator(size_t boundaryIndicator1,
			size_t boundaryIndicator2,
			shared_ptr<dealii::Triangulation<dim> > triangulation,
			dealii::Point<dim>& beginLine1, dealii::Point<dim>& endLine1,
			dealii::Point<dim>& beginLine2, dealii::Point<dim>& endLine2);


public:

	/////////////////////////////////
	// CONSTRUCTION // DESTRUCTION //
	/////////////////////////////////

	/**
	 * DEPRECATED, because boundary indicators are needed for application of boundary values
	 * @short Constructor; the periodic boundary is defined by two lines. The degrees of freedom on
	 *         the first line are forced to be equal to the opposite one on the second line, respectively.
	 *
	 *  @note  Since both lines have to match each other perfectly,
	 *         the constructor tests whether their lengths agree.
	 *
	 *  @param beginLine1 start point of line 1
	 *  @param endLine1 end point of line 1
	 *  @param beginLine2 start point of line 2
	 *  @param endLine2 end point of line 2
	 PeriodicBoundary1D(dealii::Point<2>& beginLine1, dealii::Point<2>& endLine1,
	 dealii::Point<2>& beginLine2, dealii::Point<2>& endLine2,
	 shared_ptr<dealii::Triangulation<2> > triangulation);
	 */

	/** @short Constructor; the periodic boundary is defined by two lines. The degrees of freedom on
	 *         the first line are forced to be equal to the opposite one on the second line, respectively.
	 *
	 *  @note  Since both lines have to match each other perfectly,
	 *         the constructor tests whether their lengths agree.
	 *
	 *  @param boundaryIndicator1 boundary indicator of interface line 1
	 *  @param boundaryIndicator2 boundary indicator of interface line 2
	 *  @param triangulation A (shared ptr to a) triangulation object (the mesh)
	 */
	PeriodicBoundary(size_t boundaryIndicator1, size_t boundaryIndicator2,
			shared_ptr<dealii::Triangulation<dim> > triangulation);

	/// destructor
	virtual ~PeriodicBoundary();

	/////////////////////////////////
	// APPLY BOUNDARY VALUES       //
	/////////////////////////////////

	/**
	 * @short Apply boundaries to the degrees of freedom.
	 *        This is the central function of the boundary description classes.
	 *        Periodic boundaries are put into practice by introducing constraints
	 *        x_i = x_j which force two (opposite) degrees of freedom to coincide.
	 *        In other words, degrees of freedom are eliminated.
	 *        NOTE: ELIMINATING DEGREES OF FREEDOM IS NOT POSSIBLE FOR DISCONTINUOUS GALERKIN METHODS.
	 *        When using DG methods, use PeriodicBoundary1D::getAdjacentCellAtPeriodicBoundary
	 *
	 * @note  If the discretization of the two opposite boundaries do not fit together,
	 *        a linear mapping is applied: x_i = w_j * x_j + w_k * x_k,
	 *        where x_j and x_k are the "neighboring" dofs at the opposite boundary.
	 *        NO! Linear mapping is not supported, yet!
	 *        The high-level function dofTools::make_periodicity_constraints() does not
	 *        support this.
	 *
	 * @param doFHandler The doFHandler associated with the mesh
	 * @param constraintMatrix matrix to which constraints are stored
	 */
	virtual void applyBoundaryValues(
			const shared_ptr<dealii::DoFHandler<dim> > doFHandler,
			shared_ptr<dealii::ConstraintMatrix> constraintMatrix) const;

	/**
	 * @short get the respective neighbor cell on the other side of a periodic boundary
	 *
	 * @param[in] cellID a cell ID
	 * @param[out] neighborCell the desired neighbor cell
	 *
	 * @return local face number of cell1, denoting the respective cell number
	 */
	size_t getOppositeCellAtPeriodicBoundary(const typename dealii::DoFHandler<dim>::active_cell_iterator & cell,
			typename dealii::DoFHandler<dim>::cell_iterator & neighborCell) const;

	/**
	 * @short test if a given face belongs to this boundary
	 * @param[in] cell pointer to the cell
	 * @param[in] faceBoundaryIndicator the boundary indicator of the face
	 *
	 */
	bool isFaceInBoundary(const typename dealii::DoFHandler<dim>::active_cell_iterator & cell, size_t faceBoundaryIndicator) const;

	virtual bool isPeriodic() const {
		return true;
	}

	/**
	 * @short create the map m_cells which stores the cells adjacent to the periodic boundary
	 * @short doFHandler The map is stored with doFHandler iterators in order to access degrees of freedom at the boundary.
	 */
	void createCellMap(const dealii::DoFHandler<dim>& doFHandler);

	/////////////////////////////////
	// GETTER     // SETTER        //
	/////////////////////////////////
	const dealii::Point<dim>& getBeginLine1() const {
		return m_beginLine1;
	}

	const dealii::Point<dim>& getBeginLine2() const {
		return m_beginLine2;
	}

	const dealii::Point<dim>& getEndLine1() const {
		return m_endLine1;
	}

	const dealii::Point<dim>& getEndLine2() const {
		return m_endLine2;
	}

	const shared_ptr<dealii::Triangulation<dim> >& getTriangulation() const {
		return m_triangulation;
	}

	size_t getBoundaryIndicator1() const {
		return m_boundaryIndicator1;
	}

	size_t getBoundaryIndicator2() const {
		return m_boundaryIndicator2;
	}

	const std::map<typename dealii::DoFHandler<dim>::active_cell_iterator,
			std::pair<typename dealii::DoFHandler<dim>::cell_iterator, size_t> >& getCellMap() const {
		return m_cells;
	}
};
/* PeriodicBoundary1D */

} /* namespace natrium */
#endif /* PERIODICBOUNDAR1D_H_ */
