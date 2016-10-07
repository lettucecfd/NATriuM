/*
 * GradsFunction.cpp
 *
 *  Created on: 20.09.2016
 *      Author: akraem3m
 */



#include "GradsFunction.h"

namespace natrium {

template void GradsFunction<2>(vector<double>& f, const Stencil& stencil, double rho, const dealii::Tensor<1,2>& j, const dealii::Tensor<2,2>& P);
template void GradsFunction<3>(vector<double>& f, const Stencil& stencil, double rho, const dealii::Tensor<1,3>& j, const dealii::Tensor<2,3>& P);


} /* namespace natrium */
