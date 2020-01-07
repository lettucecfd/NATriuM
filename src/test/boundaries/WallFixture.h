/*
 * WallFixture.h
 *
 *  Created on: 15.02.2018
 *      Author: akraemer
 */

#ifndef TEST_BOUNDARIES_WALLFIXTURE_H_
#define TEST_BOUNDARIES_WALLFIXTURE_H_

#include "natrium/utilities/BasicNames.h"
#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/advection/SEDGMinLee.h"
#include "natrium/advection/SemiLagrangian.h"
#include "natrium/solver/CFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"
#include "natrium/benchmarks/PeriodicTestDomain2D.h"
#include "../problemdescription/WallTestDomain2D.h"

using namespace natrium;

class BoundaryTestDensity: public dealii::Function<2> {
public:
	virtual double value(const dealii::Point<2> &, const unsigned int) const {
		return 1;
	}
};
class BoundaryTestVelocity: public dealii::Function<2> {
public:
	BoundaryTestVelocity():
			dealii::Function<2>(2){

	}
	virtual void vector_value(const dealii::Point<2> &,
			dealii::Vector<double> &values) const {
		values(0) = 0;
		values(1) = 0;
	}
};

struct WallTest {
public:
	double error_u;
	double error_v;
	WallTest(AdvectionSchemeName advection_scheme)
	{
		boost::shared_ptr<ProblemDescription<2> > problem = boost::make_shared<
				WallTestDomain2D>(1);
		boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
				SolverConfiguration>();
		configuration->setOutputDirectory("/tmp");
		configuration->setSwitchOutputOff(true);
		configuration->setNumberOfTimeSteps(100);
		configuration->setSedgOrderOfFiniteElement(2);
		configuration->setCFL(1);
		configuration->setAdvectionScheme(advection_scheme);

		CFDSolver<2> solver(configuration, problem);
		solver.run();

		// Check max velocity difference in domain
		error_u = 0.0;
		error_v = 0.0;
		double u, v = 0.0;

		const dealii::IndexSet& locally_owned_dofs (solver.getAdvectionOperator()->getLocallyOwnedDofs());
		dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
		dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());

		for (it = locally_owned_dofs.begin(); it != end; it++) {
			size_t i = *it;

			u = solver.getVelocity().at(0)(i);
			v = solver.getVelocity().at(1)(i);
			if (fabs( u - 0.01) > error_u)
				error_u = fabs(u-0.01);
			if (fabs( v - 0.0) > error_v)
				error_v = fabs(v);
		}

		error_u = dealii::Utilities::MPI::max(error_u, MPI_COMM_WORLD);
		error_v = dealii::Utilities::MPI::max(error_v, MPI_COMM_WORLD);


	}
};



#endif /* TEST_BOUNDARIES_WALLFIXTURE_H_ */
