
#ifndef LIBRARY_NATRIUM_PROBLEMDESCRIPTION_GRADSFUNCTION_H_
#define LIBRARY_NATRIUM_PROBLEMDESCRIPTION_GRADSFUNCTION_H_
/**
 * @short Definition of Grad's function to reconstruct missing distribution functions, e.g. at boundaries
 */

#include <array>

#include "deal.II/base/tensor.h"

#include "../utilities/BasicNames.h"

namespace natrium {

/**
 * @short Grad's function constructs a set of discrete distributions with respect to a certain set of moments
 * Its functional form is \f[ f_i(\rho, j, P) = w_i \left[ \rho + \frac{j_{\alpha} c_{i \alpha} }{c_s^2}
 * + \frac{1}{2 c_s^4} \left( P_{\alpha \beta} - \rho c_s^2 \delta_{\alpha \beta} \right)
 * \left( c_{i \alpha} c_{i \beta} - c_s^2 \delta_{\alpha \beta} \right) \right]. \f]
 * Checkout e.g. Dorschner et al. JCP 295 (2015) for more information.
 */
template <class Stencil, size_t dim>
void GradsFunction(vector<double>& f, const Stencil& stencil, double rho, const dealii::Tensor<1,dim>& j, const dealii::Tensor<2,dim>& P){
	const size_t Q = stencil.getQ();
	double p_part;
	double cs2 = stencil.getSpeedOfSoundSquare();
	double factor1;
	double factor2;
	assert (f.size() == Q);
	assert (rho > 0);
	for (size_t alpha = 0; alpha < Q; alpha++){
		p_part = 0;
		for (size_t i = 0; i < dim; i++){
			for (size_t j = 0; j < dim; j++){
				factor1 = P[i][j];
				factor2 = stencil.getDirection(alpha)[i]*stencil.getDirection(alpha)[j];
				if (i == j){
					factor1 -= rho * cs2;
					factor2 -= cs2;
				}
				p_part += factor1 * factor2;
			}
		}
		p_part *= (1./(2*cs2*cs2));
		f.at(alpha) = stencil.getWeight(alpha) * (rho + j * vectorToTensor<dim>(stencil.getDirection(alpha))/cs2 + p_part);
	}
}


} /* namespace natrium */

#endif /*LIBRARY_NATRIUM_PROBLEMDESCRIPTION_GRADSFUNCTION_H_*/
