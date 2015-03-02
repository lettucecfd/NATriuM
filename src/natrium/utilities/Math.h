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
inline bool is_angle_small(dealii::Point<2> vector1, dealii::Point<2> vector2,
		double thresholdDegrees = 0.0) {
	if (thresholdDegrees == 0.0) {
		return (vector1 * vector2 / (vector1.norm() * vector2.norm()) > 0.99);
	} else {
		return (vector1 * vector2 / (vector1.norm() * vector2.norm())
				> cos(thresholdDegrees / PI * 180));
	}
} /* is_angle_small */

inline double maxVelocityNorm(const vector<distributed_vector>& velocity) {
	// check sizes
	size_t dim = velocity.size();
	assert (dim > 1);
	assert (dim < 4);
	size_t n = velocity.at(0).size();
	assert (n == velocity.at(1).size());
	if (dim == 3){
		assert (n == velocity.at(2).size());
	}
	double max_norm_square = 0.0;
	for (size_t i = 0; i < n; i++) {
		double norm_square = 0.0;
		for (size_t j = 0; j < dim; j++) {
			norm_square += velocity.at(j)(i) * velocity.at(j)(i);
		}
		if (norm_square > max_norm_square) {
			max_norm_square = norm_square;
		}
	}
	return sqrt(max_norm_square);
}

inline double velocity2Norm(const vector<distributed_vector>& velocity) {
	size_t dim = velocity.size();
	assert (dim > 1);
	assert (dim < 4);
	size_t n = velocity.at(0).size();
	assert (n == velocity.at(1).size());
	if (dim == 3){
		assert (n == velocity.at(2).size());
	}

	double sum = 0.0;
	for (size_t i = 0; i < n; i++) {
		double norm_square = 0.0;
		for (size_t j = 0; j < dim; j++) {
			norm_square += velocity.at(j)(i) * velocity.at(j)(i);
		}
		sum += norm_square;
	}
	return sqrt(sum);
}


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

}
//namespace natrium

#endif /* MATH_H_ */
