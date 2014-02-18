/**
 * @file BoundaryTools.h
 * @short 
 * @date 14.11.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef BOUNDARYTOOLS_H_
#define BOUNDARYTOOLS_H_

#include <string>

#include "deal.II/base/point.h"

#include "../utilities/BasicNames.h"

namespace natrium {

namespace BoundaryTools {

/**
 * @short Check if two lines in a 2D plane are parallel and not equal to each other.
 *        If they are antiparallel, swap begin and end of the second line.
 *  @param[in] beginLine1 start point of line 1
 *  @param[in] endLine1 end point of line 1
 *  @param[in/out] beginLine2 start point of line 2
 *  @param[in/out] endLine2 end point of line 2
 *  @param[out] errorMessage contains the reason of failure when lines are not parallel
 *  @return true: if the lines are parallel or antiparallel, else: false
 */
bool checkParallelLines(const dealii::Point<2>& beginLine1,
		const dealii::Point<2>& endLine1, dealii::Point<2>& beginLine2,
		dealii::Point<2>& endLine2, std::string& errorMessage);

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
 *  @param[out] errorMessage contains the reason of failure when interfaces are no lines
 *
 *  @return true: if the interfaces are lines, false: else
 */
bool getInterfacialLinesByBoundaryIndicator(size_t boundaryIndicator1,
		size_t boundaryIndicator2,
		shared_ptr<dealii::Triangulation<2> > triangulation,
		dealii::Point<2>& beginLine1, dealii::Point<2>& endLine1,
		dealii::Point<2>& beginLine2, dealii::Point<2>& endLine2,
		std::string& errorMessage);

}

}

#endif /* BOUNDARYTOOLS_H_ */
