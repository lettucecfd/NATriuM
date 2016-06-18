/*
 * SemiLagrangianBoundaryDoFHandler.cpp
 *
 *  Created on: 10.06.2016
 *      Author: akraem3m
 */

#include "SemiLagrangianBoundaryDoFHandler.h"

namespace natrium {

/*template<size_t dim>
SemiLagrangianBoundaryDoFHandler<dim>::SemiLagrangianBoundaryDoFHandler(const Stencil& stencil) {
	// TODO Auto-generated constructor stub

}*/

template<size_t dim>
SemiLagrangianBoundaryDoFHandler<dim>::~SemiLagrangianBoundaryDoFHandler() {
	// TODO Auto-generated destructor stub
}

template class SemiLagrangianBoundaryDoFHandler<2>;
template class SemiLagrangianBoundaryDoFHandler<3>;

} /* namespace natrium */
