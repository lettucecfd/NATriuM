/*
 * ExponentialFilter_test.cpp
 *
 *  Created on: 03.02.2016
 *      Author: akraem3m
 */

#include "natrium/smoothing/ExponentialFilter.h"

#include "boost/test/unit_test.hpp"

#include "deal.II/base/polynomial_space.h"
#include "deal.II/base/polynomial.h"
#include "deal.II/fe/fe_poly.h"
#include "deal.II/fe/fe_dgp.h"
#include "deal.II/fe/fe_dgq.h"
#include "deal.II/numerics/vector_tools.h"
#include "deal.II/lac/constraint_matrix.h"
#include "deal.II/numerics/data_out.h"

#include "natrium/utilities/BasicNames.h"
#include "natrium/stencils/D2Q9.h"
#include "natrium/benchmarks/PeriodicTestDomain2D.h"
#include "natrium/advection/SEDGMinLee.h"

namespace natrium {

BOOST_AUTO_TEST_SUITE(ExponentialFilter_test)

BOOST_AUTO_TEST_CASE(ExponentialFilter_PolynomialDegree_test) {
	pout << "ExponentialFilter_GetPolynomialDegree_test..." << endl;

	// dim 1
	size_t p = 5;
	dealii::FE_DGP < 1 > fe_legendre1(p);
	dealii::Point<1> x1;
	x1(0) = 0.4;
	std::vector<dealii::Polynomials::Polynomial<double> > poly(
			dealii::Polynomials::Legendre::generate_complete_basis(p));
	for (size_t i = 0; i < p; i++) {
		BOOST_CHECK_CLOSE(fe_legendre1.shape_value(i, x1),
				poly.at(i).value(x1(0)), 1e-5);
		BOOST_CHECK_EQUAL(poly.at(i).degree(), i);
	}

	pout << "done" << endl;
}

BOOST_AUTO_TEST_CASE(ExponentialFilter_Construction2D_test) {
	pout << "ExponentialFilter_Construction2D_test..." << endl;

	size_t p = 5;
	size_t alpha = 36;
	size_t s = 2;
	size_t Nc = 1;

	dealii::QGaussLobatto<2> quadrature(p + 1);
	dealii::FE_DGQArbitraryNodes<2> fe(dealii::QGaussLobatto<1>(p + 1));

	ExponentialFilter<2> exp_filter2d(alpha, s, Nc, false, quadrature, fe);

	pout << "done" << endl;
}

BOOST_AUTO_TEST_CASE(ExponentialFilter_Construction3D_test) {
	pout << "ExponentialFilter_Construction3D_test..." << endl;

	size_t p = 5;
	size_t alpha = 36;
	size_t s = 2;
	size_t Nc = 1;

	dealii::QGaussLobatto<3> quadrature(p + 1);
	dealii::FE_DGQArbitraryNodes<3> fe(dealii::QGaussLobatto<1>(p + 1));

	ExponentialFilter<3> exp_filter3d(alpha, s, Nc, false, quadrature, fe);

	pout << "done" << endl;
}

BOOST_AUTO_TEST_CASE(ExponentialFilter_TestProjection_test) {
	pout << "ExponentialFilter_TestProjection_test..." << endl;

	size_t p = 4;
	size_t alpha = 36;
	size_t s = 2;
	size_t Nc = 1;

	dealii::QGaussLobatto<1> quadrature(p + 1);
	dealii::FE_DGQArbitraryNodes<1> fe(dealii::QGaussLobatto<1>(p + 1));

	ExponentialFilter<1> exp_filter1d(alpha, s, Nc, false, quadrature, fe);
	numeric_matrix project = exp_filter1d.getProjectToLegendre();
	numeric_vector dgq_coefficients(p + 1);
	dgq_coefficients(0) = 0.5;
	dgq_coefficients(1) = 0.6;
	dgq_coefficients(2) = 0.7;
	dgq_coefficients(3) = 0.8;
	dgq_coefficients(4) = 0.9;
	dealii::Point<1> x(0.3);
	double expected = 0.0;
	for (size_t i = 0; i < p + 1; i++) {
		expected += (dgq_coefficients(i) * fe.shape_value(i, x));
	}
	numeric_vector legendre_coefficients(p + 1);
	project.vmult(legendre_coefficients, dgq_coefficients);
	const std::vector<dealii::Polynomials::Polynomial<double> >& legendre_1d =
			exp_filter1d.getLegendre1D();
	double result = 0.0;
	for (size_t i = 0; i < p + 1; i++) {
		result += (legendre_coefficients(i) * legendre_1d.at(i).value(x(0)));
	}
	BOOST_CHECK_CLOSE(expected, result, 1e-10);

	pout << "done" << endl;
}

BOOST_AUTO_TEST_CASE(ExponentialFilter_TestFiltering_test) {
	pout << "ExponentialFilter_TestFiltering_test..." << endl;

	////////////////////////////////
	// Create a periodic function //
	////////////////////////////////
	class Modes: public dealii::Function<2> {
	private:
		size_t m_NModes;
		vector<size_t> m_frequencies;
		vector<double> m_amplitudes;
	public:
		Modes() :
				m_NModes(0) {

		}
		void addMode(size_t frequency, double amplitude) {
			m_frequencies.push_back(frequency);
			m_amplitudes.push_back(amplitude);
			m_NModes++;
		}
		virtual double value(const dealii::Point<2>& x,
				const unsigned int ) const {
			double result = 0.0;
			for (size_t i = 0; i < m_NModes; i++) {
				result += m_amplitudes.at(i)
						* (cos(2 * M_PI * m_frequencies.at(i) * x(0))
								+ sin(2 * M_PI * m_frequencies.at(i) * x(1)));
			}
			return result;
		}
	};

	Modes mode_f;
	mode_f.addMode(2, 1.0);
	mode_f.addMode(5, 1.5);
	mode_f.addMode(10, 2.0);
	mode_f.addMode(20, 0.5);

	///////////////////////////////////////
	// Project onto finite element space //
	///////////////////////////////////////
	PeriodicTestDomain2D test_domain(3);
	size_t p = 10;
	SEDGMinLee<2> sedg_operator(test_domain.getMesh(),
			test_domain.getBoundaries(), p, boost::make_shared<D2Q9>(1.0));
	sedg_operator.setupDoFs();
	sedg_operator.reassemble();
	distributed_vector vec;
	vec.reinit(sedg_operator.getSystemVector().block(0));
	dealii::VectorTools::interpolate(*sedg_operator.getDoFHandler(),
			mode_f,	vec);

	dealii::DataOut<2> data_out;
	data_out.attach_dof_handler(*sedg_operator.getDoFHandler());
	data_out.add_data_vector(vec, "original");

	// Write vtu file
	data_out.build_patches(p + 1);
	std::ofstream out_file("/tmp/natrium_original.vtu");
	data_out.write_vtu(out_file);
	out_file.close();

	///////////////////////////////////////
	// Filter and save to file //
	///////////////////////////////////////
	double alpha = 36;
	size_t s = 2;
	size_t Nc = 1;
	ExponentialFilter<2> exp_filter(alpha, s, Nc, false, *sedg_operator.getQuadrature(),
			*sedg_operator.getFe());
	exp_filter.applyFilter(*sedg_operator.getDoFHandler(), vec);
	std::ofstream out_file2("/tmp/natrium_filtered.vtu");
	dealii::DataOut<2> data_out2;
	data_out2.attach_dof_handler(*sedg_operator.getDoFHandler());
	data_out2.add_data_vector(vec, "filtered");
	data_out2.build_patches( p + 1 );
	data_out2.write_vtu(out_file2);
	out_file2.close();

	pout << "done" << endl;
}

BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
