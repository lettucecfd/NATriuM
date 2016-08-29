/*
 * sl-stability.cpp
 *
 *  Created on: 24.08.2016
 *      Author: akraem3m
 */

#include <fstream>

#include "boost/shared_ptr.hpp"

#include "deal.II/fe/fe_dgq.h"
#include "deal.II/lac/lapack_full_matrix.h"

#include "natrium/utilities/BasicNames.h"

using std::cout;
using std::endl;
using std::ofstream;

// one-dimensional sl advection on a single cell (periodic)
double spectral_radius(boost::shared_ptr<dealii::FiniteElement<1, 1> > fe,
		double dx) {

	assert(fe->has_support_points());
	assert(dx <= 1);
	const std::vector<dealii::Point<1> > support_points(
			fe->get_unit_support_points());
	size_t n_points = support_points.size();
	assert(n_points > 0);
	assert(n_points == fe->n_dofs_per_cell());

	// build spatial operator in 1d
	dealii::LAPACKFullMatrix<double> m(n_points, n_points);
	m = 0;
	dealii::Point<1> p;
	for (size_t i = 0; i < n_points; i++) {
		p(0) = support_points.at(i)(0) - dx;
		// periodic boundary
		if (p(0) < 0)
			p(0) += 1;
		for (size_t j = 0; j < n_points; j++) {
			m(i, j) += fe->shape_value(j, p);
		}
	}

	// compute spectrum
	m.compute_eigenvalues();
	// write eigenvalues to vector
	double abs_max_eigenvalue = 0.0;
	for (size_t i = 0; i < m.n_cols(); i++) {
		if (abs(m.eigenvalue(i)) > abs_max_eigenvalue) {
			abs_max_eigenvalue = abs(m.eigenvalue(i));
		}
	}
	return abs_max_eigenvalue;

}

int main() {

	// FEs
	std::vector<boost::shared_ptr<dealii::FiniteElement<1, 1> > > fes;
	for (size_t p = 1; p < 11; p++) {
		// FE_DGQ GLL
		fes.push_back(
				boost::make_shared<dealii::FE_DGQArbitraryNodes<1> >(
						dealii::QGaussLobatto<1>(p + 1)));
	}
	for (size_t p = 1; p < 11; p++) {
		// FE_DGQ
		fes.push_back(boost::make_shared<dealii::FE_DGQ<1> >(p));
	}
	for (size_t p = 1; p < 11; p++) {
		// FE_DGQ
		fes.push_back(
				boost::make_shared<dealii::FE_DGQArbitraryNodes<1> >(
						dealii::QGaussChebyshev<1>(p + 1)));
	}
	for (size_t p = 1; p < 11; p++) {
		// FE_DGQ
		fes.push_back(
				boost::make_shared<dealii::FE_DGQArbitraryNodes<1> >(
						dealii::QGaussLobattoChebyshev<1>(p + 1)));
	}
	for (size_t p = 1; p < 11; p++) {
		// FE_DGQ
		fes.push_back(
				boost::make_shared<dealii::FE_DGQArbitraryNodes<1> >(
						dealii::QGaussRadauChebyshev<1>(p + 1)));
	}


	// make legend file
	std::stringstream f1;
	f1 << getenv("NATRIUM_HOME") << "/semilagrange-stability/legend.txt";
	std::stringstream f2;
	f2 << getenv("NATRIUM_HOME")
			<< "/semilagrange-stability/spectral_radius.txt";

	ofstream legend(f1.str());
	ofstream results(f2.str());

	for (size_t i = 0; i < fes.size(); i++) {
		legend << i + 1 << " " << fes.at(i)->get_name() << endl;
	}
	legend.close();

	// compute stability
	const size_t n_steps = 1000;
	for (size_t j = 0; j <= n_steps; j++) {
		double dx = (1.0 * j) / n_steps;
		results << dx;
		for (size_t i = 0; i < fes.size(); i++) {
			results << " " << spectral_radius(fes.at(i), dx);
		}
		results << endl;
	}
	results.close();

}
