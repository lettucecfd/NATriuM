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
#include "deal.II/grid/grid_tools.h"
#include "deal.II/dofs/dof_handler.h"
#include "deal.II/lac/dynamic_sparsity_pattern.h"//compressed_sparsity_pattern.h"

#include "Boundary.h"
#include "../utilities/NATriuMException.h"
#include "../utilities/Logging.h"
#include "../utilities/DealiiExtensions.h"


namespace natrium {

/**
 * @short Exception for impossible periodic boundary,
 * e.g. when the interfaces don't have the same length.
 */
class PeriodicBoundaryNotPossible: public NATriuMException {
private:
	std::string message;
public:
	PeriodicBoundaryNotPossible(const char *msg,
			const std::stringstream & additionalInfo = std::stringstream()) :
			NATriuMException(msg), message(msg) {
		LOG(DETAILED) << "Additional information on error: "
				<< additionalInfo.str().c_str() << endl;

	}
	PeriodicBoundaryNotPossible(const string& msg,
			const std::stringstream & additionalInfo = std::stringstream()) :
			NATriuMException(msg), message(msg) {
		LOG(DETAILED) << "Additional information on error: "
				<< additionalInfo.str().c_str() << endl;
	}
	~PeriodicBoundaryNotPossible() throw () {
	}
	const char *what() const throw () {
		return this->message.c_str();
	}
};

/**
 * @short  A periodic boundary condition. Periodic boundaries have to be set before the grid is refined!
 * @note   First use in step-1 tutorial.
 */
template<size_t dim>
class PeriodicBoundary: public Boundary<dim> {
	friend class own_double_less_periodic;
private:

	/// triangulation object
	boost::shared_ptr<Mesh<dim> > m_triangulation;

	/// boundary indicator of first interfacial line
	size_t m_boundaryIndicator1;

	/// boundary indicator of second interfacial line
	size_t m_boundaryIndicator2;

	/// direction of the boundary
	size_t m_direction;

	/// Container for all cells that belong to this boundary
	/// stored as <accessor to cell, (accessor to opposite cell, boundary face at opposite cell) > /
	PeriodicCellMap<dim> m_cells;

	dealii::DoFHandler<dim> const * m_doFHandler;

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
	 boost::shared_ptr<Mesh<2> > triangulation);
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
			size_t direction, boost::shared_ptr<Mesh<dim> > triangulation);

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
	size_t getOppositeCellAtPeriodicBoundary(
			const typename dealii::DoFHandler<dim>::active_cell_iterator & cell,
			typename dealii::DoFHandler<dim>::cell_iterator & neighborCell) const;

	/**
	 * @short test if a given face belongs to this boundary
	 * @param[in] cell pointer to the cell
	 * @param[in] face_id local face id (< 6)
	 *
	 */
	bool isFaceInBoundary(
			const typename dealii::DoFHandler<dim>::active_cell_iterator & cell,
			size_t face_id) const;


	virtual bool isPeriodic() const {
		return true;
	}
	virtual bool isLinearFluxBoundary() const {
		return false;
	}
	virtual bool isSLBoundary() const {
		return false;
	}

	/**
	 * @short create the map m_cells which stores the cells adjacent to the periodic boundary
	 * @param doFHandler The map is stored with doFHandler iterators in order to access degrees of freedom at the boundary.
	 */
	void createCellMap(const dealii::DoFHandler<dim>& doFHandler);

	/** @short
	 * check if all cells in the cell map have appropriate boundary indicators
	 */
	void checkCellMap();

	/**
	 * @short This function does nothing; just to satisfy the interface.
	 * 		  The Periodic Boundary conditions are directly incorporated in make_sparser_flux_sparsity_pattern
	 */
	void addToSparsityPattern(dealii::BlockDynamicSparsityPattern&, size_t,
			size_t, size_t) const;

	/**
	 * Transform a point into its equivalent across the periodic boundary
	 */
	dealii::Point<dim> coordinatesAcrossPeriodicBoundary(
			const dealii::Point<dim>& p,
			const typename dealii::DoFHandler<dim>::active_cell_iterator& cell);
	//TODO should be inlined




/////////////////////////////////
// GETTER     // SETTER        //
/////////////////////////////////

	const boost::shared_ptr<Mesh<dim> >& getMesh() const {
		return m_triangulation;
	}

	size_t getBoundaryIndicator1() const {
		return m_boundaryIndicator1;
	}

	size_t getBoundaryIndicator2() const {
		return m_boundaryIndicator2;
	}

	const PeriodicCellMap<dim>& getCellMap() const {
		return m_cells;
	}

	size_t getDirection() const {
		return m_direction;
	}

	const boost::shared_ptr<Mesh<dim> >& getTriangulation() const {
		return m_triangulation;
	}


}
;
/* PeriodicBoundary1D */

} /* namespace natrium */
#endif /* PERIODICBOUNDAR1D_H_ */
