/**
 * @file PeriodicBoundary1D.h
 * @short Description of a periodic boundary on a line.
 * @date 25.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef PERIODICBOUNDARY1D_H_
#define PERIODICBOUNDARY1D_H_

#include "deal.II/base/point.h"

#include "BoundaryDescription.h"

namespace natrium {

/**
 * @short  A periodic boundary condition on a line, appropriate for 2D problem descriptions.
 * @note   First use in step-1 tutorial.
 */
class PeriodicBoundary1D: public BoundaryDescription<1> {
private:

	/// start point of line 1
	dealii::Point<2> m_beginLine1;

	/// end point of line 1
	dealii::Point<2> m_endLine1;

	/// start point of line 2
	dealii::Point<2> m_beginLine2;

	/// end point of line 2
	dealii::Point<2> m_endLine2;

public:

	/** @short Constructor; the periodic boundary is defined by two lines. The degrees of freedom on
	 *         the first line are forced to be equal to the opposite one on the second line, respectively.
	 *
	 *  @note  Since both lines have to match each other perfectly,
	 *         the constructor tests whether their lengths agree.
	 *
	 *  @param beginLine1 start point of line 1
	 *  @param endLine1 end point of line 1
	 *  @param beginLine2 start point of line 2
	 *  @param endLine2 end point of line 2
	 */
	PeriodicBoundary1D(dealii::Point<2> beginLine1, dealii::Point<2> endLine1,
			dealii::Point<2> beginLine2, dealii::Point<2> endLine2);

	/// destructor
	virtual ~PeriodicBoundary1D();

	/**
	 * @short Apply boundaries to the degrees of freedom.
	 *        This is the central function of the boundary description classes.
	 *        Periodic boundaries are put into practice by introducing constraints
	 *        x_i = x_j which force two (opposite) degrees of freedom to coincide.
	 *
	 * @note  The user has to make shure that the discretization of two opposite lines are identical
	 *        (position of vertices). This is because no mapping between nodes is available at the
	 *        current development stage of the code.
	 *        If positions do not coincide applyBoundaryValues throws an error.
	 *
	 * @param triangulation A triangulation object (the mesh)
	 * @param doFHandler The doFHandler associated with the mesh
	 */
	virtual void applyBoundaryValues(
			shared_ptr<dealii::Triangulation<2> > triangulation,
			shared_ptr<dealii::DoFHandler<2> > doFHandler);

};
/* PeriodicBoundary1D */

} /* namespace natrium */
#endif /* PERIODICBOUNDAR1D_H_ */
