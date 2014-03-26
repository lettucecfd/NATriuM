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
#include <map>

#include "deal.II/base/point.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"
#include "deal.II/dofs/dof_handler.h"
#include <deal.II/lac/compressed_sparsity_pattern.h>

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
	 * @param doFHandler The map is stored with doFHandler iterators in order to access degrees of freedom at the boundary.
	 */
	void createCellMap(const dealii::DoFHandler<dim>& doFHandler);

	/**
	 * @short modify sparsity pattern so that the fluxes over periodic boundary can be incorporated
	 * @param cSparse the block-sparsity pattern
	 * @param n_blocks the number of blocks in cSparse
	 * @param n_dofs_per_row number of degrees of freedom per block (normally: overall degrees of freedom on grid)
	 * @param dofs_per_cell number of degrees of freedom per cell
	 */
	void addToSparsityPattern(dealii::BlockCompressedSparsityPattern&  cSparse, size_t n_blocks, size_t n_dofs_per_block, size_t dofs_per_cell) const;


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
