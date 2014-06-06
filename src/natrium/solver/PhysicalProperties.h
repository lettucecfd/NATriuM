/**
 * @file SolverConfiguration.h
 * @short Calculation of physical properties
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef PHYSICALPROPERTIES_H_
#define PHYSICALPROPERTIES_H_

#include "../utilities/BasicNames.h"

namespace natrium {

template<size_t dim> class PhysicalProperties {
public:

	/// constructor
	PhysicalProperties();

	/// destructor
	virtual ~PhysicalProperties();

	// kinetic energy
	static double kineticEnergy(const vector<distributed_vector>& u, const distributed_vector& rho);

};

} /* namespace natrium */

#endif /* PHYSICALPROPERTIES_H_ */
