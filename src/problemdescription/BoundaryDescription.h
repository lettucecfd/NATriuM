/**
 * @file BoundaryDescription.h
 * @short Abstract class for Description of a Boundary object
 * @date 24.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef BOUNDARYDESCRIPTION_H_
#define BOUNDARYDESCRIPTION_H_

# include "../utilities/BasicNames.h"

namespace natrium {

template<size_t dim> class BoundaryDescription {

public:
	BoundaryDescription();
	virtual ~BoundaryDescription();
};


template<size_t dim>
inline natrium::BoundaryDescription<dim>::BoundaryDescription() {
}

template<size_t dim>
inline natrium::BoundaryDescription<dim>::~BoundaryDescription() {
}

} /* namespace natrium */

#endif /* BOUNDARYDESCRIPTION_H_ */
