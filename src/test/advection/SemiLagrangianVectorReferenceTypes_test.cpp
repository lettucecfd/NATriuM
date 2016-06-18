/*
 * SemiLagrangianVectorReferenceTypes.cpp
 *
 *  Created on: 10.06.2016
 *      Author: akraem3m
 */

#include "natrium/advection/SemiLagrangianVectorReferenceTypes.h"

#include "deal.II/dofs/dof_handler.h"
#include "deal.II/fe/fe_dgq.h"
#include "deal.II/base/quadrature.h"
#include "natrium/advection/SemiLagrangianBoundaryDoFHandler.h"
#include "natrium/advection/SemiLagrangian.h"
#include "natrium/problemdescription/BoundaryCollection.h"
#include "natrium/benchmarks/PeriodicTestDomain2D.h"

#include "boost/test/unit_test.hpp"

namespace natrium {

BOOST_AUTO_TEST_SUITE(SemiLagrangianVectorReferenceTypes_test)

BOOST_AUTO_TEST_CASE(GeneralizedDoF_Construction_test) {
	pout << "GeneralizedDoF_Construction_test" << endl;

	GeneralizedDoF d(true, 1, 1);
	BOOST_CHECK_EQUAL(d.getAlpha(), size_t(1));
	BOOST_CHECK_EQUAL(d.getIndex(), size_t(1));
	BOOST_CHECK_EQUAL(d.isSecondaryBoundaryHit(), true);

	GeneralizedDoF e(d);
	BOOST_CHECK_EQUAL(d.getAlpha(), e.getAlpha());
	BOOST_CHECK_EQUAL(d.getIndex(), e.getIndex());
	BOOST_CHECK_EQUAL(d.isSecondaryBoundaryHit(), e.isSecondaryBoundaryHit());

	pout << "done." << endl;
} /* GeneralizedDoF_Construction_test*/

BOOST_AUTO_TEST_CASE(SecondaryBoundaryDoFVector_Construction_test) {
	pout << "SecondaryBoundaryDoFVector_Construction_test" << endl;

	size_t fe_order = 1;
	size_t refinementLevel = 3;
	PeriodicTestDomain2D periodic(refinementLevel);
	periodic.refineAndTransform();
	dealii::DoFHandler<2> dof_handler(*periodic.getMesh());
	dealii::FE_DGQArbitraryNodes<2> fe(
			dealii::QGaussLobatto<1>(fe_order + 1));
	dof_handler.distribute_dofs(fe);

	SecondaryBoundaryDoFVector sbdv(dof_handler.locally_owned_dofs());

	dof_handler.clear();
	pout << "done." << endl;
} /* SecondaryBoundaryDoFVector_Construction_test*/

BOOST_AUTO_TEST_CASE(SecondaryBoundaryDoFVector_appendBoundaryDoFAndCompress_test) {
	pout << "SecondaryBoundaryDoFVector_appendBoundaryDoFAndCompress_test" << endl;

	size_t fe_order = 1;
	size_t refinementLevel = 3;
	PeriodicTestDomain2D periodic(refinementLevel);
	periodic.refineAndTransform();
	dealii::DoFHandler<2> dof_handler(*periodic.getMesh());
	dealii::FE_DGQArbitraryNodes<2> fe(
			dealii::QGaussLobatto<1>(fe_order + 1));
	dof_handler.distribute_dofs(fe);

	SecondaryBoundaryDoFVector sbdv(dof_handler.locally_owned_dofs());

	GeneralizedDoF a = sbdv.appendSecondaryBoundaryDoF();
	BOOST_CHECK_EQUAL(a.getIndex(),
			dof_handler.locally_owned_dofs().nth_index_in_set(0));
	BOOST_CHECK_EQUAL(a.getAlpha(), size_t(0));
	BOOST_CHECK_EQUAL(a.isSecondaryBoundaryHit(), true);
	BOOST_CHECK_EQUAL(sbdv.getSecondaryBoundaryIndices().n_elements(), size_t(1));

	size_t n = 2;
	for (size_t i = 0; i < n; i++) {
		sbdv.appendSecondaryBoundaryDoF();
	}

	BOOST_CHECK_EQUAL(sbdv.getSecondaryBoundaryIndices().n_elements(), n + 1);

	sbdv.compress();

	BOOST_CHECK_EQUAL(sbdv.getSecondaryBoundaryValues().size(), n+1);

	dof_handler.clear();
	pout << "done." << endl;
} /* SecondaryBoundaryDoFVector_appendBoundaryDoFAndCompress_test*/

BOOST_AUTO_TEST_CASE(SemiLagrangianVectorAccess_operatorSubscript_test) {
	pout << "SemiLagrangianVectorAccess_operatorSubscript_test" << endl;

	size_t fe_order = 1;
	size_t refinementLevel = 3;
	PeriodicTestDomain2D periodic(refinementLevel);
	periodic.refineAndTransform();
	dealii::DoFHandler<2> dof_handler(*periodic.getMesh());
	dealii::FE_DGQArbitraryNodes<2> fe(
			dealii::QGaussLobatto<1>(fe_order + 1));
	dof_handler.distribute_dofs(fe);

	vector<distributed_vector> vec;
	for (size_t i = 0; i < 9; i++) {
		distributed_vector tmp(dof_handler.locally_owned_dofs(),
				MPI_COMM_WORLD);
		vec.push_back(tmp);
	}
	DistributionFunctions f_old (vec);

	vector<distributed_vector> vec2;
	for (size_t i = 0; i < 9; i++) {
		distributed_vector tmp(dof_handler.locally_owned_dofs(),
				MPI_COMM_WORLD);
		vec2.push_back(tmp);
	}
	DistributionFunctions f_new (vec2);

	SecondaryBoundaryDoFVector sbdv(dof_handler.locally_owned_dofs());

	SemiLagrangianVectorAccess f(f_old, f_new, sbdv);

	GeneralizedDoF a = sbdv.appendSecondaryBoundaryDoF();
	GeneralizedDoF b = sbdv.appendSecondaryBoundaryDoF();
	GeneralizedDoF c(false, 1, 1);
	GeneralizedDoF d(true, 0, 0);

	sbdv.compress();

	f[a] = 2.0;
	f[b] = 3.14;
	f[c] = 2.7;

	BOOST_CHECK_EQUAL(f[a], 2.0);
	BOOST_CHECK_EQUAL(f[b], 3.14);
	BOOST_CHECK_EQUAL(f[c], 2.7);

	BOOST_CHECK_EQUAL(f(a), 2.0);
	BOOST_CHECK_EQUAL(f(b), 3.14);
	BOOST_CHECK_EQUAL(f(c), 2.7);


	dof_handler.clear();

	pout << "done." << endl;
} /* SemiLagrangianVectorAccess_operatorSubscript_test*/

BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
