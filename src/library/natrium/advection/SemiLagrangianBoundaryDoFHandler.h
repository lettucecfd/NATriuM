/*
 * SemiLagrangianBoundaryDoFHandler.h
 *
 *  Created on: 10.06.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_ADVECTION_SEMILAGRANGIANBOUNDARYDOFHANDLER_H_
#define LIBRARY_NATRIUM_ADVECTION_SEMILAGRANGIANBOUNDARYDOFHANDLER_H_

#include "deal.II/base/index_set.h"
#include "deal.II/base/tensor.h"
#include "deal.II/dofs/dof_handler.h"

#include "BoundaryHit.h"
#include "../problemdescription/BoundaryCollection.h"
#include "../utilities/BasicNames.h"


namespace natrium {

/**
 * @short
 */
template<size_t dim>
class SemiLagrangianBoundaryDoFHandler {
private:

	dealii::IndexSet m_dofsThatDependOnBoundaryValues;

	distributed_vector m_boundaryValues;

	std::vector<BoundaryHit<dim> > m_boundaryPoints;

	distributed_vector m_secondaryBoundaryValues;

	std::vector<BoundaryHit<dim> > m_secondaryBoundaryPoints;

	// boost::shared_ptr<BoundaryCollection<dim> >& m_boundaries;

public:
	SemiLagrangianBoundaryDoFHandler();
	virtual ~SemiLagrangianBoundaryDoFHandler();

	/**
	 * @short
	 */
	void addBoundaryHit(BoundaryHit<dim>& boundary_hit){

	}
	void calculateBoundaryValues(double time_of_next_step){}
	void applyBoundaryValues(distributed_block_vector& f){}



};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_ADVECTION_SEMILAGRANGIANBOUNDARYDOFHANDLER_H_ */
