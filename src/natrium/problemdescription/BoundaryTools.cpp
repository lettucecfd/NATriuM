/**
 * @file BoundaryTools.cpp
 * @short 
 * @date 14.11.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "BoundaryTools.h"

bool natrium::BoundaryTools::checkParallelLines(
		const dealii::Point<2>& beginLine1, const dealii::Point<2>& endLine1,
		dealii::Point<2>& beginLine2, dealii::Point<2>& endLine2,
		std::string& errorMessage) {

	// check input
	double lengthLine1 = beginLine1.distance(endLine1);
	double lengthLine2 = beginLine2.distance(endLine2);

	// assert that points are not equal to one another
	if ((lengthLine1 == 0) or (lengthLine2 == 0)
			or (beginLine1.distance(beginLine2) == 0)
			or (endLine1.distance(endLine2) == 0)
			or (beginLine1.distance(endLine2) == 0)
			or (beginLine2.distance(endLine1) == 0)) {
		errorMessage.clear();
		errorMessage.append(
				"Two of the points defining a periodic boundary are equal. That is not allowed.");
		return false;
	}

	// assert that both lines have same length (up to 1%)
	if (abs(lengthLine1 - lengthLine2) / lengthLine1 > 0.01) {
		errorMessage.clear();
		errorMessage.append(
				"The two lines defining a periodic boundary must have the same length.");
		return false;
	}

	// assert that interfaces are parallel (anything else would need different handling)
	dealii::Point<2> differenceVector1 = endLine1 - beginLine1;
	dealii::Point<2> differenceVector2 = endLine2 - beginLine2;

	if (not Math::is_angle_small(differenceVector1, differenceVector2)) {
		// try to fix the problem by swapping begin and end
		dealii::Point<2> tmp = beginLine2;
		differenceVector2 = beginLine2 - endLine2;
		if (not Math::is_angle_small(differenceVector1, differenceVector2)) {
			// return
			errorMessage.clear();
			errorMessage.append(
					"The two lines defining a periodic boundary must be parallel to each other.");
			return false;
		}
		// re-orient line 2
		beginLine2 = endLine2;
		endLine2 = tmp;
	}
	return true;
}
