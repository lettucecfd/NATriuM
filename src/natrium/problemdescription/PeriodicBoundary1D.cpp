/**
 * @file PeriodicBoundary1D.cpp
 * @short Description of a periodic boundary on a line.
 * @date 25.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "PeriodicBoundary1D.h"

namespace natrium {



PeriodicBoundary1D::PeriodicBoundary1D(dealii::Point<2> beginLine1,
		dealii::Point<2> endLine1, dealii::Point<2> beginLine2,
		dealii::Point<2> endLine2) {
}

PeriodicBoundary1D::~PeriodicBoundary1D() {
}

void PeriodicBoundary1D::applyBoundaryValues(
		shared_ptr<dealii::Triangulation<2> > triangulation,
		shared_ptr<dealii::DoFHandler<2> > doFHandler) {
}

} /* namespace natrium */
