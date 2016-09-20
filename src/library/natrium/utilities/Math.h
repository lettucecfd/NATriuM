/**
 * @file Math.h
 * @short Definition of basic math functionality
 * @date 18.02.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef MATH_H_
#define MATH_H_

#include <functional>
#include <math.h>

#include "deal.II/base/index_set.h"
#include "deal.II/base/tensor.h"

#include "BasicNames.h"

namespace natrium {

/// namespace which contains basic math functions
namespace Math {

const double PI = 3.141592653589793238462;

/// scalar product
inline double scalar_product(const numeric_vector& x, const numeric_vector& y) {
	return x * y;
}

/// scale existing vector
inline void scale_vector(double a, numeric_vector& x) {
	x *= a;
}

/// scalar times vector
inline numeric_vector scalar_vector(double a, const numeric_vector& x) {
	numeric_vector y(x);
	y *= a;
	return y;
}

// add vectors
inline void add_vector(numeric_vector& x, const numeric_vector& y) {
	x += y;
}

// subtract vectors
inline void subtract_vector(numeric_vector& x, const numeric_vector& y) {
	x -= y;
}

// 2-norm
inline double euclidean_norm(numeric_vector& x) {
	return x.l2_norm();
}

// check if the angle between to 2d vectors is small
// The formula stems from the cosine theorem: <a,b> / (|a| |b|) = cos( angle(a,b) ).
// Asserting abs[ <a,b> / (|a| |b|) ] > 0.99 is equivalent to:    angle(a,b) < 8 degrees.
inline bool is_angle_small(dealii::Tensor<1, 2> vector1,
		dealii::Tensor<1, 2> vector2, double thresholdDegrees = 0.0) {
	if (thresholdDegrees == 0.0) {
		return (vector1 * vector2 / (vector1.norm() * vector2.norm()) > 0.99);
	} else {
		return (vector1 * vector2 / (vector1.norm() * vector2.norm())
				> cos(thresholdDegrees / PI * 180));
	}
} /* is_angle_small */

/**
 * @short Calculate the maximum of the euclidean velocity norm
 * @param[in] velocity Velocity vector.
 * @param[in] locally_owned_dofs Index set that contains the locally owned degrees of freedom.
 * @note not implemented for block vectors, because the velocity should not come as a block vector in natrium
 */
template<class VECTOR>
inline double maxVelocityNorm(const vector<VECTOR>& velocity,
		const dealii::IndexSet& locally_owned_dofs) {
	// check sizes
	size_t dim = velocity.size();
	assert(dim > 1);
	assert(dim < 4);
	size_t n = velocity.at(0).size();
	assert(n == velocity.at(1).size());
	if (dim == 3) {
		assert(n == velocity.at(2).size());
	}
	double max_norm_square = 0.0;
	//for all degrees of freedom on current processor
	dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
	dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());
	for (; it != end; it++) {
		size_t i = *it;
		double norm_square = 0.0;
		for (size_t j = 0; j < dim; j++) {
			norm_square += velocity.at(j)(i) * velocity.at(j)(i);
		}
		if (norm_square > max_norm_square) {
			max_norm_square = norm_square;
		}
	}
	double global_mpi = dealii::Utilities::MPI::min_max_avg(max_norm_square,
			MPI_COMM_WORLD).max;
	return sqrt(global_mpi);
}

/**
 * @short Calculate the euclidean norm of the velocity over all points (not in the finite element space, but element-wise)
 * @param[in] velocity Velocity vector.
 * @param[in] locally_owned_dofs Index set that contains the locally owned degrees of freedom.
 * @note not implemented for block vectors, because the velocity should not come as a block vector in natrium
 * @note In contrast to maxVelocityNorm, this function is seperately implemented for distributed vectors,
 *       (specialized template instantiation) because a global sum would corrupt the solution, when using local vectors.
 *       This does not allow to declare the function as inline.
 */
template< class VECTOR>
double velocity2Norm(const vector<VECTOR>& velocity,
		const dealii::IndexSet& locally_owned_dofs);

} /* Math */

/**
 * @short function to compare doubles as map keys;
 *        regards two doubles equal if they are in a epsilon(=1e-7)-range.
 */
class own_double_less: public std::binary_function<double, double, bool> {
public:
	own_double_less(double arg_ = 1e-6) :
			epsilon(arg_) {
	}
	bool operator()(const double &left, const double &right) const {
		// you can choose other way to make decision
		// (The original version is: return left < right;)
		return (fabs(left - right) > epsilon) && (left < right);
	}
	double epsilon;
};

template<size_t dim>
dealii::Tensor<1, dim> vectorToTensor(const numeric_vector& v) {
	assert(v.size() == dim);
	dealii::Tensor<1, dim> t;
	t[0] = v(0);
	t[1] = v(1);
	if (dim == 3)
		t[2] = v(2);
	return t;
}

}
//namespace natrium

#endif /* MATH_H_ */
