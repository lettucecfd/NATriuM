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
		size_t boundaryIndicator,
		shared_ptr<dealii::Function<dim> > boundaryDensity,
		shared_ptr<dealii::Function<dim> > boundaryVelocity) :
		m_boundaryIndicator(boundaryIndicator), m_boundaryDensity(
				boundaryDensity), m_boundaryVelocity(boundaryVelocity) {
}
template MinLeeBoundary<2>::MinLeeBoundary(size_t boundaryIndicator,
		shared_ptr<dealii::Function<2> > boundaryDensity,
		shared_ptr<dealii::Function<2> > boundaryVelocity);
template MinLeeBoundary<3>::MinLeeBoundary(size_t boundaryIndicator,
		shared_ptr<dealii::Function<3> > boundaryDensity,
		shared_ptr<dealii::Function<3> > boundaryVelocity);

template<size_t dim>
MinLeeBoundary<dim>::MinLeeBoundary(size_t boundaryIndicator,
		const dealii::Vector<double>& velocity) :
		m_boundaryIndicator(boundaryIndicator), m_boundaryDensity(
				make_shared<BoundaryDensity>()), m_boundaryVelocity(
				make_shared<BoundaryVelocity>(velocity)) {
}
template MinLeeBoundary<2>::MinLeeBoundary(size_t boundaryIndicator,
		const dealii::Vector<double>& velocity);
/*template MinLeeBoundary<3>::MinLeeBoundary(size_t boundaryIndicator,
 const dealii::Vector<double>& velocity);
 */
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
			continue;

		// add
		for (size_t I = 0; I < n_blocks; I++) {
			for (size_t J = 0; J < n_blocks; J++) {
				if (I
						== boltzmannModel.getIndexOfOppositeDirection(J + 1)
								- 1) {
					// get global degrees of freedom
					cell->get_dof_indices(localDoFIndices);
					for (size_t i = 0; i < n_dofs_per_cell; i++) {
						for (size_t j = 0; j < n_dofs_per_cell; j++) {
							// TODO: Efficient implementation (for SEDG, MinLee Boundary): only at diagonal
							cSparse.block(I, J).add(localDoFIndices.at(i),
									localDoFIndices.at(j));
							//cSparse.block(J, I).add(localDoFIndices.at(i),
							//		localDoFIndices.at(j));
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

template<size_t dim> void MinLeeBoundary<dim>::assembleBoundary(size_t alpha,
		const typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
		size_t faceNumber, dealii::FEFaceValues<dim>& feFaceValues,
		const BoltzmannModel& boltzmannModel,
		const std::map<size_t, size_t>& q_index_to_facedof,
		distributed_sparse_block_matrix& systemMatrix,
		distributed_block_vector& systemVector, bool useCentralFlux) const {
	// let the feFaceValues object calculate all the values needed at the boundary
	feFaceValues.reinit(cell, faceNumber);
	const vector<double> &JxW = feFaceValues.get_JxW_values();
	const vector<dealii::Point<dim> > &normals =
			feFaceValues.get_normal_vectors();

// TODO clean up construction; or (better): assembly directly
	dealii::FullMatrix<double> cellFaceMatrix(feFaceValues.dofs_per_cell);
	cellFaceMatrix = 0;
	dealii::FullMatrix<double> oppositeDirectionCellFaceMatrix(
			feFaceValues.dofs_per_cell);
	oppositeDirectionCellFaceMatrix = 0;
	dealii::Vector<double> cellFaceVector(feFaceValues.dofs_per_cell);
	cellFaceVector = 0;

	// calculate prefactor
	for (size_t q = 0; q < feFaceValues.n_quadrature_points; q++) {
		size_t thisDoF = q_index_to_facedof.at(q);

		// get density and velocity at boundary point
		double density = m_boundaryDensity->value(
				feFaceValues.quadrature_point(q));
		dealii::Vector<double> velocity(dim);
		m_boundaryVelocity->vector_value(feFaceValues.quadrature_point(q),
				velocity);

		// calculate matrix entries
		double prefactor = JxW.at(q);
		double exn = 0.0;
		// calculate scalar product
		for (size_t i = 0; i < dim; i++) {		// TODO efficient multiplication
			exn += normals.at(q)(i) * boltzmannModel.getDirection(alpha)(i);
		}
		prefactor *= exn;

		// calculate vector entry
		double exu = 0.0;
		// calculate scalar product
		for (size_t i = 0; i < dim; i++) {		// TODO efficient multiplication
			exu += velocity(i) * boltzmannModel.getDirection(alpha)(i);
		}

		if (useCentralFlux) {
			cellFaceMatrix(thisDoF, thisDoF) = 0.5 * prefactor;
			oppositeDirectionCellFaceMatrix(thisDoF, thisDoF) = 0.5 * prefactor;
			cellFaceVector(thisDoF) = -prefactor
					* boltzmannModel.getWeight(alpha) * density * exu
					/ boltzmannModel.getSpeedOfSoundSquare();
		} else if (exn < 0) {
			cellFaceMatrix(thisDoF, thisDoF) = prefactor;
			oppositeDirectionCellFaceMatrix(thisDoF, thisDoF) = -prefactor;
			cellFaceVector(thisDoF) = -2 * prefactor
					* boltzmannModel.getWeight(alpha) * density * exu
					/ boltzmannModel.getSpeedOfSoundSquare();
		}
	}

// get DoF indices
// TODO cut out construction (allocation); Allocating two vectors in most inner loop is too expensive
	vector<dealii::types::global_dof_index> localDoFIndices(
			feFaceValues.dofs_per_cell);
	cell->get_dof_indices(localDoFIndices);

/// Distribute to global matrix
/// TODO Loop only over diagonal
	for (size_t i = 0; i < feFaceValues.dofs_per_cell; i++) {
		for (size_t j = 0; j < feFaceValues.dofs_per_cell; j++) {
			systemMatrix.block(alpha - 1, alpha - 1).add(localDoFIndices[i],
					localDoFIndices[j], cellFaceMatrix(i, j));
			systemMatrix.block(alpha - 1,
					boltzmannModel.getIndexOfOppositeDirection(alpha) - 1).add(
					localDoFIndices[i], localDoFIndices[j],
					oppositeDirectionCellFaceMatrix(i, j));
		}
		systemVector.block(alpha - 1)(localDoFIndices[i]) += cellFaceVector(i);
	}
}
template void MinLeeBoundary<2>::assembleBoundary(size_t alpha,
		const typename dealii::DoFHandler<2>::active_cell_iterator& cell,
		size_t faceNumber, dealii::FEFaceValues<2>& feFaceValues,
		const BoltzmannModel& boltzmannModel,
		const std::map<size_t, size_t>& q_index_to_facedof,
		distributed_sparse_block_matrix& systemMatrix,
		distributed_block_vector& systemVector, bool useCentralFlux) const;
template void MinLeeBoundary<3>::assembleBoundary(size_t alpha,
		const typename dealii::DoFHandler<3>::active_cell_iterator& cell,
		size_t faceNumber, dealii::FEFaceValues<3>& feFaceValues,
		const BoltzmannModel& boltzmannModel,
		const std::map<size_t, size_t>& q_index_to_facedof,
		distributed_sparse_block_matrix& systemMatrix,
		distributed_block_vector& systemVector, bool useCentralFlux) const;

} /* namespace natrium */
