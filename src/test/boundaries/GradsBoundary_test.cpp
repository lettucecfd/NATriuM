#include "natrium/boundaries/GradsBoundary.h"

#include <mpi.h>
#include <string>

#include "boost/test/unit_test.hpp"

#include "deal.II/base/tensor.h"
#include "deal.II/grid/grid_generator.h"
#include "deal.II/base/function.h"
#include "deal.II/base/function_lib.h"
#include "deal.II/dofs/dof_handler.h"
#include "deal.II/fe/fe_dgq.h"
#include "deal.II/base/quadrature_lib.h"


#include "natrium/stencils/D2Q9.h"
#include "natrium/stencils/D3Q15.h"
#include "natrium/stencils/D3Q19.h"
#include "natrium/stencils/D3Q27.h"
#include "natrium/solver/DistributionFunctions.h"
#include "natrium/solver/CFDSolver.h"
#include "natrium/benchmarks/CouetteFlowGrad2D.h"
#include "natrium/boundaries/SLBoundary.h"

#include "natrium/problemdescription/BoundaryCollection.h"
#include "natrium/advection/SemiLagrangian.h"

#include "natrium/utilities/BasicNames.h"

using namespace natrium;

BOOST_AUTO_TEST_SUITE(GradsBoundary_test)

BOOST_AUTO_TEST_CASE(GradsBoundary_Constructor_test) {

	pout << "GradsBoundary_Constructor_test..." << endl;

	boost::shared_ptr<dealii::Function<2> > f = boost::make_shared<dealii::ConstantFunction<2> >(0.1, 2);
	boost::shared_ptr<dealii::Function<3> > f3d = boost::make_shared<dealii::ConstantFunction<3> >(0.1, 3);

	GradsBoundary<2> boundary(0,f);
	GradsBoundary<3> boundary2(0,f3d);

	//GradsBoundary<2,PRESCRIBED_PRESSURE> boundary3(0,f);
	//GradsBoundary<3,PRESCRIBED_PRESSURE> boundary4(0,f3d);

	pout << "done." << endl;

} /* GradsFunction_Construction_test */

BOOST_AUTO_TEST_CASE(GradsBoundary_Velocity2D_test) {

	pout << "GradsBoundary_Velocity2D_test..." << endl;

	boost::shared_ptr<dealii::Function<2> > f1 = boost::make_shared<dealii::ConstantFunction<2> >(0.1, 2);
	CouetteFlowGrad2D couette(0.1, 0.01, 3);
	couette.refineAndTransform();

	// distribute dofs
	SemiLagrangian<2> sl(couette,2,boost::make_shared<D2Q9>(),0.001);
	sl.setupDoFs();

	DistributionFunctions f;
	distributed_vector rho;
	vector<distributed_vector> u;
	f.reinit(9, sl.getDoFHandler()->locally_owned_dofs(), MPI_COMM_WORLD);
	f.compress(dealii::VectorOperation::add);
	rho.reinit(sl.getDoFHandler()->locally_owned_dofs(), MPI_COMM_WORLD);
	rho.compress(dealii::VectorOperation::add);
	for (size_t i = 0; i < 2; i++){
		distributed_vector tmp;
		tmp.reinit(sl.getDoFHandler()->locally_owned_dofs(), MPI_COMM_WORLD);
		tmp.compress(dealii::VectorOperation::add);
		u.push_back(tmp);
	}
	assert (u.size() == 2);
	assert (f.size() > 0);
	assert (f.at(0).size() > 0);
	assert (rho.size() > 0);
	assert (u.at(0).size() > 0);


	double beta = 1;
	D2Q9 d2q9;
	//boundary->apply(f, rho, u, sl, beta, d2q9);

	// TODO: Calculate u beforehand

	// TODO: Test afterwards

	pout << "done." << endl;

} /* GradsBoundary_Velocity2D_test */


BOOST_AUTO_TEST_CASE(GradsBoundary_Couette2D_test) {
	pout << "GradsBoundary_Couette2D_test..." << endl;

	double viscosity = 1e-1;
	double u0 = 0.1;
	size_t ref_level = 3;
	boost::shared_ptr<ProblemDescription<2> > couette = boost::make_shared<CouetteFlowGrad2D>(viscosity, u0, ref_level);

	// TODO: finalize
	pout << "done." << endl;
} /* GradsBoundary_Velocity2D_test */


BOOST_AUTO_TEST_SUITE_END()

