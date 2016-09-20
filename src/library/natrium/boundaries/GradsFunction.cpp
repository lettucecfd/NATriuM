/*
 * GradsFunction.cpp
 *
 *  Created on: 20.09.2016
 *      Author: akraem3m
 */



#include "GradsFunction.h"

namespace natrium {

class D2Q9;
class D3Q15;
class D3Q19;
class D3Q27;

template void GradsFunction<D2Q9, 2>(vector<double>& f, const D2Q9& stencil, double rho, const dealii::Tensor<1,2>& j, const dealii::Tensor<2,2>& P);
template void GradsFunction<D3Q15, 3>(vector<double>& f, const D3Q15& stencil, double rho, const dealii::Tensor<1,3>& j, const dealii::Tensor<2,3>& P);
template void GradsFunction<D3Q19, 3>(vector<double>& f, const D3Q19& stencil, double rho, const dealii::Tensor<1,3>& j, const dealii::Tensor<2,3>& P);
template void GradsFunction<D3Q27, 3>(vector<double>& f, const D3Q27& stencil, double rho, const dealii::Tensor<1,3>& j, const dealii::Tensor<2,3>& P);


} /* namespace natrium */
