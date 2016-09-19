/*
 * LinearBoundaryRhoU.cpp
 *
 *  Created on: 08.12.2015
 *      Author: akraem3m
 */

#include "LinearFluxBoundaryRhoU.h"

namespace natrium {

template<size_t dim>
LinearFluxBoundaryRhoU<dim>::LinearFluxBoundaryRhoU(size_t boundaryIndicator,
		boost::shared_ptr<dealii::Function<dim> > boundaryDensity,
		boost::shared_ptr<dealii::Function<dim> > boundaryVelocity) :
		LinearFluxBoundary<dim>(boundaryIndicator, boundaryDensity,
				boundaryVelocity,
				BoundaryTools::COUPLE_ONLY_OPPOSITE_DISTRIBUTIONS,
				BoundaryTools::COUPLE_ONLY_SINGLE_POINTS) {

}

/// constructor
template<size_t dim>
LinearFluxBoundaryRhoU<dim>::LinearFluxBoundaryRhoU(size_t boundaryIndicator,
		const dealii::Vector<double>& velocity) :
		LinearFluxBoundaryRhoU(boundaryIndicator,
				boost::make_shared<BoundaryTools::BoundaryDensity<dim> >(),
				boost::make_shared<BoundaryTools::BoundaryVelocity<dim> >(velocity)) {

}

template<size_t dim>
LinearFluxBoundaryRhoU<dim>::~LinearFluxBoundaryRhoU() {
}

template<size_t dim> void LinearFluxBoundaryRhoU<dim>::assembleBoundary(
		size_t alpha,
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
	const vector<dealii::Tensor<1, dim> > &normals =
			feFaceValues.get_all_normal_vectors();

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
		double density = LinearFluxBoundary<dim>::getBoundaryDensity()->value(
				feFaceValues.quadrature_point(q));
		dealii::Vector<double> velocity(dim);
		LinearFluxBoundary<dim>::getBoundaryVelocity()->vector_value(
				feFaceValues.quadrature_point(q), velocity);

		// calculate matrix entries
		double prefactor = JxW.at(q);
		double exn = 0.0;
		// calculate scalar product
		for (size_t i = 0; i < dim; i++) {	// TODO efficient multiplication
			exn += normals.at(q)[i] * stencil.getDirection(alpha)(i);
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

// Explicit instantiation
template class LinearFluxBoundaryRhoU<2> ;
template class LinearFluxBoundaryRhoU<3> ;

} /* namespace natrium */
