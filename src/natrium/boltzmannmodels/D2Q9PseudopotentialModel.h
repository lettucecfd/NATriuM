/*
 * D2Q9PseudopotentialModel.h
 *
 *  Created on: Nov 11, 2014
 *      Author: kk
 */

#ifndef D2Q9PSEUDOPOTENTIALMODEL_H_
#define D2Q9PSEUDOPOTENTIALMODEL_H_

#include "D2Q9Model.h"

#include "deal.II/grid/tria.h"
#include "deal.II/numerics/vector_tools.h"

#include "../utilities/BasicNames.h"

#include "AdvectionOperator.h"

namespace natrium {

	class D2Q9PseudopotentialModel: public D2Q9Model {
	public:

		/// constructor
		D2Q9PseudopotentialModel(double scaling,
				shared_ptr<SolverConfiguration> configuration,
				shared_ptr<Triangulation<dim>> triangulation,
				DistributionFunctions& f,
				distributed_vector& densities);

		/// destructor
		virtual ~D2Q9PseudopotentialModel();

		virtual double getEquilibriumDistribution(size_t i, const numeric_vector& u, const double rho = 1) const;

		virtual vector<double> getDensityGradient() ;

		virtual void getInteractionForce(shared_ptr<SolverConfiguration> configuration, size_t i, const numeric_vector& u, const double rho = 1) ;
	};

} /* namespace natrium */

#endif /* D2Q9PSEUDOPOTENTIALMODEL_H_ */
