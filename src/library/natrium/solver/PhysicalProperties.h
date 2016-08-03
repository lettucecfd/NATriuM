/**
 * @file SolverConfiguration.h
 * @short Calculation of physical properties
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef PHYSICALPROPERTIES_H_
#define PHYSICALPROPERTIES_H_

#include "../advection/AdvectionOperator.h"
#include "../utilities/BasicNames.h"

namespace natrium {

template<size_t dim> class PhysicalProperties {
public:

	/// constructor
	PhysicalProperties();

	/// destructor
	virtual ~PhysicalProperties();

	// kinetic energy
	static double kineticEnergy(const vector<distributed_vector>& u, const distributed_vector& rho,
			boost::shared_ptr<AdvectionOperator<dim> > advection);

	static double enstrophy(const vector<distributed_vector>& u,
			boost::shared_ptr<AdvectionOperator<dim> > advection);

	static double mass(const distributed_vector& rho, boost::shared_ptr<AdvectionOperator<dim> > advection);


	/// Pressure
	static double maximalPressure(const distributed_vector& rho, const double speedOfSound, double & minimalPressure);

	/// Flow factor
	static double meanVelocityX(const distributed_vector& ux, boost::shared_ptr<AdvectionOperator<dim> > advection);

	static double entropy(const DistributionFunctions& f, boost::shared_ptr<AdvectionOperator<dim> > advection);

};

} /* namespace natrium */


#endif /* PHYSICALPROPERTIES_H_ */
