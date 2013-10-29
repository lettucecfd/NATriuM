/**
 * @file PeriodicBoundary1D.h
 * @short Description of a periodic boundary on a line.
 * @date 25.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef PERIODICBOUNDARY1D_H_
#define PERIODICBOUNDARY1D_H_

#include "exception"
#include "string"

#include "deal.II/base/point.h"

#include "BoundaryDescription.h"

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
 * @short  A periodic boundary condition on a line, appropriate for 2D problem descriptions.
 * @note   First use in step-1 tutorial.
 */
class PeriodicBoundary1D: public BoundaryDescription<1> {
private:

	/// boundary indicator of first interfacial line
	size_t m_boundaryIndicator1;

	/// boundary indicator of second interfacial line
	size_t m_boundaryIndicator2;

	/// start point of line 1
	dealii::Point<2> m_beginLine1;

	/// end point of line 1
	dealii::Point<2> m_endLine1;

	/// start point of line 2
	dealii::Point<2> m_beginLine2;

	/// end point of line 2
	dealii::Point<2> m_endLine2;

	/// triangulation object
	shared_ptr<dealii::Triangulation<2> > m_triangulation;

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
			shared_ptr<dealii::Triangulation<2> > triangulation,
			dealii::Point<2>& beginLine1, dealii::Point<2>& endLine1,
			dealii::Point<2>& beginLine2, dealii::Point<2>& endLine2);

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
	PeriodicBoundary1D(size_t boundaryIndicator1, size_t boundaryIndicator2,
			shared_ptr<dealii::Triangulation<2> > triangulation);

	/// destructor
	virtual ~PeriodicBoundary1D();

	/////////////////////////////////
	// APPLY BOUNDARY VALUES       //
	/////////////////////////////////

	/**
	 * @short Apply boundaries to the degrees of freedom.
	 *        This is the central function of the boundary description classes.
	 *        Periodic boundaries are put into practice by introducing constraints
	 *        x_i = x_j which force two (opposite) degrees of freedom to coincide.
	 *
	 * @note  If the discretization of the two opposite boundaries do not fit together,
	 *        a linear mapping is applied: x_i = w_j * x_j + w_k * x_k,
	 *        where x_j and x_k are the "neighboring" dofs at the opposite boundary.
	 *
	 * @param doFHandler The doFHandler associated with the mesh
	 * @param constraintMatrix matrix to which constraints are stored
	 */
	virtual void applyBoundaryValues(
			const shared_ptr<dealii::DoFHandler<2> > doFHandler,
			shared_ptr<dealii::ConstraintMatrix> constraintMatrix) const;

	/////////////////////////////////
	// GETTER     // SETTER        //
	/////////////////////////////////

	const dealii::Point<2>& getBeginLine1() const {
		return m_beginLine1;
	}

	const dealii::Point<2>& getBeginLine2() const {
		return m_beginLine2;
	}

	const dealii::Point<2>& getEndLine1() const {
		return m_endLine1;
	}

	const dealii::Point<2>& getEndLine2() const {
		return m_endLine2;
	}

	const shared_ptr<dealii::Triangulation<2> >& getTriangulation() const {
		return m_triangulation;
	}

};
/* PeriodicBoundary1D */

} /* namespace natrium */
#endif /* PERIODICBOUNDAR1D_H_ */
