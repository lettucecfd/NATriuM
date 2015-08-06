/*
 * MinLeeBoundary.cpp
 *
 *  Created on: 26.03.2014
 *      Author: kraemer
 */

#include "MinLeeBoundary.h"

#include "deal.II/lac/compressed_sparsity_pattern.h"
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>

namespace natrium {

template class BoundaryDensity<2> ;
template class BoundaryDensity<3> ;
template class BoundaryVelocity<2> ;
template class BoundaryVelocity<3> ;

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
				make_shared<BoundaryDensity<dim> >()), m_boundaryVelocity(
				make_shared<BoundaryVelocity<dim> >(velocity)) {
}
template MinLeeBoundary<2>::MinLeeBoundary(size_t boundaryIndicator,
		const dealii::Vector<double>& velocity);
template MinLeeBoundary<3>::MinLeeBoundary(size_t boundaryIndicator,
		const dealii::Vector<double>& velocity);

template<size_t dim> void MinLeeBoundary<dim>::addToSparsityPattern(
#ifdef WITH_TRILINOS_MPI
		dealii::TrilinosWrappers::SparsityPattern& cSparse,
#else
		dealii::DynamicSparsityPattern& cSparse,
#endif
		const dealii::DoFHandler<dim>& doFHandler,
		const Stencil& stencil) const {

	// ConstraintMatrix can be used for a more efficient distribution to global sparsity patterns
	const dealii::ConstraintMatrix constraints;

	size_t n_dofs_per_cell = doFHandler.get_fe().dofs_per_cell;
	std::vector<dealii::types::global_dof_index> dofs_on_this_cell(
			n_dofs_per_cell);
	std::vector<dealii::types::global_dof_index> dofs_on_this_face;
	dofs_on_this_face.reserve(n_dofs_per_cell);

	// couple opposite distribution functions at boundaries
	// iterate over all cells
	typename dealii::DoFHandler<dim>::active_cell_iterator cell =
			doFHandler.begin_active();
	typename dealii::DoFHandler<dim>::active_cell_iterator endc =
			doFHandler.end();
	for (; cell != endc; ++cell) {
		if (cell->is_locally_owned()) {
			for (size_t i = 0; i < dealii::GeometryInfo<dim>::faces_per_cell;
					i++) {
				if (cell->face(i)->at_boundary()) {
					if (cell->face(i)->boundary_indicator()
							== m_boundaryIndicator) {
						cell->get_dof_indices(dofs_on_this_cell);
						dofs_on_this_face.clear();
						for (size_t j = 0; j < n_dofs_per_cell; j++) {
							if (cell->get_fe().has_support_on_face(j, i)) {
								dofs_on_this_face.push_back(
										dofs_on_this_cell.at(j));
							}
						}

						// add
						// couple only individual dofs with each other
						std::vector<dealii::types::global_dof_index> dof_this(
								1);
						std::vector<dealii::types::global_dof_index> dof_this2(
								1);
						for (size_t i = 0; i < dofs_on_this_face.size(); i++) {
							dof_this.at(0) = dofs_on_this_face.at(i);
							dof_this2.at(0) = dofs_on_this_face.at(i);
							// Add entries to sparsity pattern
							constraints.add_entries_local_to_global(dof_this,
									dof_this2, cSparse, true);
						}
					} /* end if boundary indicator */
				} /* end if at boundary */
			} /* end for all faces */
		} /* end if is locally owned */
	} /* end forall cells */
}

template void MinLeeBoundary<2>::addToSparsityPattern(
#ifdef WITH_TRILINOS_MPI
		dealii::TrilinosWrappers::SparsityPattern& cSparse,
#else
		dealii::DynamicSparsityPattern& cSparse,
#endif
		const dealii::DoFHandler<2>& doFHandler, const Stencil& stencil) const;
template void MinLeeBoundary<3>::addToSparsityPattern(
#ifdef WITH_TRILINOS_MPI
		dealii::TrilinosWrappers::SparsityPattern& cSparse,
#else
		dealii::DynamicSparsityPattern& cSparse,
#endif
		const dealii::DoFHandler<3>& doFHandler, const Stencil& stencil) const;


template<size_t dim> void MinLeeBoundary<dim>::assembleBoundary(size_t alpha,
		const typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
		size_t faceNumber, dealii::FEFaceValues<dim>& feFaceValues,
		const Stencil& stencil,
		const std::map<size_t, size_t>& q_index_to_facedof,
		const vector<double> & inverseLocalMassMatrix,
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
		for (size_t i = 0; i < dim; i++) {	// TODO efficient multiplication
			exn += normals.at(q)(i) * stencil.getDirection(alpha)(i);
		}
		prefactor *= exn;

		// calculate vector entry
		double exu = 0.0;
		// calculate scalar product
		for (size_t i = 0; i < dim; i++) {	// TODO efficient multiplication
			exu += velocity(i) * stencil.getDirection(alpha)(i);
		}

		if (useCentralFlux) {
			cellFaceMatrix(thisDoF, thisDoF) = 0.5 * prefactor;
			oppositeDirectionCellFaceMatrix(thisDoF, thisDoF) = 0.5 * prefactor;
			cellFaceVector(thisDoF) = -prefactor * stencil.getWeight(alpha)
					* density * exu / stencil.getSpeedOfSoundSquare();
		} else if (exn < 0) {
			cellFaceMatrix(thisDoF, thisDoF) = prefactor;
			oppositeDirectionCellFaceMatrix(thisDoF, thisDoF) = -prefactor;
			cellFaceVector(thisDoF) = -2 * prefactor * stencil.getWeight(alpha)
					* density * exu / stencil.getSpeedOfSoundSquare();
		}
	}

// get DoF indices
// TODO cut out construction (allocation); Allocating two vectors in most inner loop is too expensive
	vector<dealii::types::global_dof_index> localDoFIndices(
			feFaceValues.dofs_per_cell);
	cell->get_dof_indices(localDoFIndices);

/// Distribute to global matrix
	for (size_t i = 0; i < feFaceValues.dofs_per_cell; i++) {
		if (cell->get_fe().has_support_on_face(i, faceNumber)) {
			systemMatrix.block(alpha - 1, alpha - 1).add(localDoFIndices[i],
					localDoFIndices[i],
					cellFaceMatrix(i, i) * inverseLocalMassMatrix.at(i));
			systemMatrix.block(alpha - 1,
					stencil.getIndexOfOppositeDirection(alpha) - 1).add(
					localDoFIndices[i], localDoFIndices[i],
					oppositeDirectionCellFaceMatrix(i, i)
							* inverseLocalMassMatrix.at(i));
			systemVector.block(alpha - 1)(localDoFIndices[i]) +=
					(cellFaceVector(i) * inverseLocalMassMatrix.at(i));
		}
	}
}
template void MinLeeBoundary<2>::assembleBoundary(size_t alpha,
		const typename dealii::DoFHandler<2>::active_cell_iterator& cell,
		size_t faceNumber, dealii::FEFaceValues<2>& feFaceValues,
		const Stencil& stencil,
		const std::map<size_t, size_t>& q_index_to_facedof,
		const vector<double> & inverseLocalMassMatrix,
		distributed_sparse_block_matrix& systemMatrix,
		distributed_block_vector& systemVector, bool useCentralFlux) const;
template void MinLeeBoundary<3>::assembleBoundary(size_t alpha,
		const typename dealii::DoFHandler<3>::active_cell_iterator& cell,
		size_t faceNumber, dealii::FEFaceValues<3>& feFaceValues,
		const Stencil& stencil,
		const std::map<size_t, size_t>& q_index_to_facedof,
		const vector<double> & inverseLocalMassMatrix,
		distributed_sparse_block_matrix& systemMatrix,
		distributed_block_vector& systemVector, bool useCentralFlux) const;

} /* namespace natrium */
