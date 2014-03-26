/*
 * MinLeeBoundary.cpp
 *
 *  Created on: 26.03.2014
 *      Author: kraemer
 */

#include <problemdescription/MinLeeBoundary.h>

#include <deal.II/dofs/dof_handler.h>

namespace natrium {

template<size_t dim> MinLeeBoundary<dim>::MinLeeBoundary(
		size_t boundaryIndicator) :
		m_boundaryIndicator(boundaryIndicator) {
}
template MinLeeBoundary<2>::MinLeeBoundary(size_t boundaryIndicator);
template MinLeeBoundary<3>::MinLeeBoundary(size_t boundaryIndicator);

template<size_t dim> void MinLeeBoundary<dim>::addToSparsityPattern(
		dealii::BlockCompressedSparsityPattern& cSparse,
		const dealii::DoFHandler<dim>& doFHandler,
		const BoltzmannModel& boltzmannModel) const {
	size_t n_blocks = boltzmannModel.getQ() - 1;
	size_t n_dofs_per_cell = doFHandler.get_fe().dofs_per_cell;
	std::vector<dealii::types::global_dof_index> localDoFIndices(
			n_dofs_per_cell);

	// couple opposite distribution functions at boundaries
	// iterate over all cells
	typename dealii::DoFHandler<dim>::active_cell_iterator cell =
			doFHandler.begin_active();
	typename dealii::DoFHandler<dim>::active_cell_iterator endc =
			doFHandler.end();
	for (; cell != endc; ++cell) {
		// check if at least one face is at this wall
		bool isAtThisWall = false;
		for (size_t i = 0; i < dealii::GeometryInfo<2>::faces_per_cell; i++) {
			if (cell->face(i)->at_boundary()) {
				if (cell->face(i)->boundary_indicator()
						== m_boundaryIndicator) {
					isAtThisWall = true;
					break;
				}
			}
		}
		if (not isAtThisWall)
			break;

		// add
		for (size_t I = 0; I < n_blocks; I++) {
			for (size_t J = I + 1; J < n_blocks; J++) {
				if (I == boltzmannModel.getIndexOfOppositeDirection(J)) {
					// get global degrees of freedom
					cell->get_dof_indices(localDoFIndices);
					for (size_t i = 0; i < n_dofs_per_cell; i++) {
						for (size_t j = 0; j < n_dofs_per_cell; j++) {
							// Efficient implementation (for SEDG, MinLee Boundary): only at diagonal
							cSparse.block(I, J).add(localDoFIndices.at(i),
									localDoFIndices.at(j));
							cSparse.block(J, I).add(localDoFIndices.at(i),
									localDoFIndices.at(j));
						}
					}
				} /* end if opposite boundary*/
			} /* end block J */
		} /* end block I */
	} /* end forall cells */
}
template void MinLeeBoundary<2>::addToSparsityPattern(
		dealii::BlockCompressedSparsityPattern& cSparse,
		const dealii::DoFHandler<2>& doFHandler,
		const BoltzmannModel& boltzmannModel) const;
template void MinLeeBoundary<3>::addToSparsityPattern(
		dealii::BlockCompressedSparsityPattern& cSparse,
		const dealii::DoFHandler<3>& doFHandler,
		const BoltzmannModel& boltzmannModel) const;

} /* namespace natrium */
