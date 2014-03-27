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

template<size_t dim> void MinLeeBoundary<dim>::assembleBoundary(
		size_t alpha,
		const typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
		size_t faceNumber, dealii::FEFaceValues<dim>& feFaceValues,
		const BoltzmannModel& boltzmannModel,
		const std::map<size_t, size_t>& q_index_to_facedof,
		distributed_sparse_block_matrix& systemMatrix,
		distributed_block_vector& systemVector) const {
	// let the feFaceValues object calculate all the values needed at the boundary
	feFaceValues.reinit(cell, faceNumber);
	const vector<double> &JxW = feFaceValues.get_JxW_values();
	const vector<dealii::Point<dim> > &normals = feFaceValues.get_normal_vectors();

// TODO clean up construction; or (better): assembly directly
	dealii::FullMatrix<double> cellFaceMatrix(feFaceValues.dofs_per_cell);
	cellFaceMatrix = 0;
	dealii::Vector<double> cellFaceVector(feFaceValues.dofs_per_cell);

	// calculate prefactor
	for (size_t q = 0; q < feFaceValues.n_quadrature_points; q++) {
			size_t thisDoF = q_index_to_facedof.at(q);

			// calculate matrix entries
			vector<double> factor(JxW);
			double exn = 0.0;
			// calculate scalar product
			for (size_t i = 0; i < dim; i++) {		// TODO efficient multiplication
				exn += normals.at(q)(i)
						* boltzmannModel.getDirection(alpha)(i);
			}
			for (size_t i = 0; i < factor.size(); i++) {
				factor.at(i) *= exn;
			}

			//
			double density = 1/0;
			cellFaceMatrix(thisDoF, thisDoF) = factor.at(q);
			//cellFaceVector(thisDoF) = -2* factor.at(q) * boltzmannModel.getWeight(direction) * density ;
	}

	// get DoF indices
	// TODO cut out construction (allocation); Allocating two vectors in most inner loop is too expensive
		vector<dealii::types::global_dof_index> localDoFIndices(
				feFaceValues.dofs_per_cell);
		cell->get_dof_indices(localDoFIndices);

	/// Distribute to global matrix
		for (size_t i = 0; i < feFaceValues.dofs_per_cell; i++) {
			for (size_t j = 0; j < feFaceValues.dofs_per_cell; j++) {
				systemMatrix.block(alpha-1, alpha-1).add(localDoFIndices[i],
						localDoFIndices[j], cellFaceMatrix(i, j));
				systemMatrix.block(alpha, boltzmannModel.getIndexOfOppositeDirection(alpha)-1).add(localDoFIndices[i],
						localDoFIndices[j], -cellFaceMatrix(i, j));
			}
		}


}
template void MinLeeBoundary<2>::assembleBoundary(size_t alpha,
		const typename dealii::DoFHandler<2>::active_cell_iterator& cell,
		size_t faceNumber, dealii::FEFaceValues<2>& feFaceValues,
		const BoltzmannModel& boltzmannModel,
		const std::map<size_t, size_t>& q_index_to_facedof,
		distributed_sparse_block_matrix& systemMatrix,
		distributed_block_vector& systemVector) const;
template void MinLeeBoundary<3>::assembleBoundary(size_t alpha,
		const typename dealii::DoFHandler<3>::active_cell_iterator& cell,
		size_t faceNumber, dealii::FEFaceValues<3>& feFaceValues,
		const BoltzmannModel& boltzmannModel,
		const std::map<size_t, size_t>& q_index_to_facedof,
		distributed_sparse_block_matrix& systemMatrix,
		distributed_block_vector& systemVector) const;

} /* namespace natrium */
