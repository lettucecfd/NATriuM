/*
 * DealiiExtensions.h
 *
 *  Created on: Dec 11, 2014
 *      Author: kraemer
 */

#ifndef DEALIIEXTENSIONS_H_
#define DEALIIEXTENSIONS_H_


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/table.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/dofs/function_map.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/hp/mapping_collection.h>

#include <vector>
#include <set>
#include <map>



namespace natrium {

namespace DealIIExtensions{

//DEAL_II_NAMESPACE_OPEN

using namespace dealii;


/**
 * @short Like dealii::DoFTools::make_flux_sparsity_pattern but does only create
 * 		  non-zero entries for the DoFs situated on faces
 */
template <class DH, class SparsityPattern>
  void
  make_sparser_flux_sparsity_pattern (const DH                  &dof,
                              SparsityPattern           &sparsity,
                              const ConstraintMatrix    &constraints,
                              const bool                keep_constrained_dofs = true,
                              const types::subdomain_id  subdomain_id = numbers::invalid_unsigned_int);



template <class DH, class SparsityPattern>
 void
 make_sparser_flux_sparsity_pattern (const DH        &dof,
                             SparsityPattern &sparsity);


//DEAL_II_NAMESPACE_CLOSE

} /* namespace DealIIExtensions */

} /* namespace natrium */



#endif /* DEALIIEXTENSIONS_H_ */


