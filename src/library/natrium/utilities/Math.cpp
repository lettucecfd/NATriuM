/*
 * Math.cpp
 *
 *  Created on: 16.11.2015
 *      Author: akraem3m
 */

#include "Math.h"

namespace natrium {
namespace Math {

template
double maxVelocityNorm<numeric_vector>(const vector<numeric_vector>& velocity,
		const dealii::IndexSet& locally_owned_dofs);
template<>
double velocity2Norm<numeric_vector>(const vector<numeric_vector>& velocity,
		const dealii::IndexSet& locally_owned_dofs) {
	size_t dim = velocity.size();
	assert(dim > 1);
	assert(dim < 4);
	size_t n = velocity.at(0).size();
	assert(n == velocity.at(1).size());
	if (dim == 3) {
		assert(n == velocity.at(2).size());
	}

	double sum = 0.0;
	//for all degrees of freedom on current processor
	dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
	dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());
	for (; it != end; it++) {
		size_t i = *it;
		double norm_square = 0.0;
		for (size_t j = 0; j < dim; j++) {
			norm_square += velocity.at(j)(i) * velocity.at(j)(i);
		}
		sum += norm_square;
	}
	return sqrt(sum);
}

template
double maxVelocityNorm<distributed_vector>(const vector<distributed_vector>& velocity,
		const dealii::IndexSet& locally_owned_dofs);
template<>
double velocity2Norm<distributed_vector>(const vector<distributed_vector>& velocity,
		const dealii::IndexSet& locally_owned_dofs) {
	size_t dim = velocity.size();
	assert(dim > 1);
	assert(dim < 4);
	size_t n = velocity.at(0).size();
	assert(n == velocity.at(1).size());
	if (dim == 3) {
		assert(n == velocity.at(2).size());
	}

	double sum = 0.0;
	//for all degrees of freedom on current processor
	dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
	dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());
	for (; it != end; it++) {
		size_t i = *it;
		double norm_square = 0.0;
		for (size_t j = 0; j < dim; j++) {
			norm_square += velocity.at(j)(i) * velocity.at(j)(i);
		}
		sum += norm_square;
	}
	// the only difference to the specialized template instatiation for numeric_vector
	double global_mpi =
			dealii::Utilities::MPI::min_max_avg(sum, MPI_COMM_WORLD).sum;
	return sqrt(global_mpi);
}



} /* namespace Math */
} /* namespace natrium */
