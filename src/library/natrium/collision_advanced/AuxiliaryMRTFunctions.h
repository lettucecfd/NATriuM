#ifndef LIBRARY_NATRIUM_COLLISION_AUXILIARY_MRT_FUNCTIONS_H_
#define LIBRARY_NATRIUM_COLLISION_AUXILIARY_MRT_FUNCTIONS_H_

#include "../stencils/Stencil.h"
#include "../utilities/BasicNames.h"
#include "../utilities/ConfigNames.h"
#include <array>

using std::array;

namespace natrium {

namespace AuxiliaryMRTFunctions {

/**
 * @short Exception class for Collision
 */
class MRTException: public NATriuMException {
private:
	std::string message;
public:
	MRTException(const char *msg) :
			NATriuMException(msg), message(msg) {
	}
	MRTException(const string& msg) :
			NATriuMException(msg), message(msg) {
	}
	~MRTException() throw () {
	}
	const char *what() const throw () {
		return this->message.c_str();
	}
};

// M: moment transform
// T: inverse moment transform
// diag: matrix of relaxation times

// Weighted MRT basis by Dellar (2003)
class MRTDellarD2Q9 {
	static const array<array<double, 9>, 9> moment_trafo;
	static const array<array<double, 9>, 9> inverse_trafo;
};

// Cartesian MRT basis by Lallemand and Luo (2000)
class MRTLallemandD2Q9 {
	static const array<array<double, 9>, 9> moment_trafo;
	static const array<array<double, 9>, 9> inverse_trafo;
};

// Cartesian MRT basis by D'Humieres et al. (2002)
class MRTDHumieresD3Q19 {
	static const array<array<double, 19>, 19> moment_trafo;
	static const array<array<double, 19>, 19> inverse_trafo;
};

// factory functions
template<size_t T_Q>
const array<array<double, T_Q>, T_Q> make_M(MomentBasis basis);
template<size_t T_Q>
const array<array<double, T_Q>, T_Q> make_T(MomentBasis basis);
template<size_t T_Q>
const array<double, T_Q> make_diag(double tau, MomentBasis basis,
		RelaxMode relax_mode);

// helper functions
template<size_t T_m, size_t T_n>
void matrix_vector_product(const array<array<double, T_n>, T_m> matrix,
		const array<double, T_n>& source, array<double, T_m>& destination) {
	destination = {};
	for (size_t i = 0; i < T_m; i++) {
		for (size_t j = 0; j < T_n; j++) {
			destination[i] += matrix[i][j] * source[j];
		}
	}
}

template<size_t T_m, size_t T_n, size_t T_k>
void matrix_matrix_product(const array<array<double, T_n>, T_m> A,
		const array<array<double, T_k>, T_n> B,
		array<array<double, T_k>, T_m> result) {
	result = {};
	for (size_t i = 0; i < T_m; i++) {
		for (size_t j = 0; j < T_n; j++) {
			for (size_t k = 0; k < T_k; k++) {
				result[i][k] += A[i][j] * B[j][k];
			}
		}
	}
}

}
/* namespace AuxiliaryMRTFunctions */

} /* namespace natrium */

#endif /* include guard */
