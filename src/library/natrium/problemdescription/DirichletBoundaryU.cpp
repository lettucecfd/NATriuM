/*
 * DirichletBoundaryU.cpp
 *
 *  Created on: 08.12.2015
 *      Author: akraem3m
 */

#include "DirichletBoundaryU.h"

namespace natrium {

template<size_t dim>
DirichletBoundaryU<dim>::DirichletBoundaryU(size_t boundaryIndicator,
		boost::shared_ptr<dealii::Function<dim> > boundaryVelocity) :
		DirichletBoundary<dim>(boundaryIndicator, NULL,
				boundaryVelocity,
				BoundaryTools::COUPLE_ALL_DISTRIBUTIONS,
				BoundaryTools::COUPLE_ONLY_SINGLE_POINTS) {
}

/// constructor
template<size_t dim>
DirichletBoundaryU<dim>::DirichletBoundaryU(size_t boundaryIndicator,
		const dealii::Vector<double>& velocity):
		DirichletBoundaryU<dim>(boundaryIndicator,
				boost::make_shared<BoundaryVelocity<dim> >(velocity)){
}

template<size_t dim>
DirichletBoundaryU<dim>::~DirichletBoundaryU() {
}

template<size_t dim> void DirichletBoundaryU<dim>::assembleBoundary(
		size_t alpha,
		const typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
		size_t faceNumber, dealii::FEFaceValues<dim>& feFaceValues,
		const Stencil& stencil,
		const std::map<size_t, size_t>& q_index_to_facedof,
		const vector<double> & inverseLocalMassMatrix,
		distributed_sparse_block_matrix& systemMatrix,
		distributed_block_vector& , bool useCentralFlux) const {
	// let the feFaceValues object calculate all the values needed at the boundary
	feFaceValues.reinit(cell, faceNumber);
	const vector<double> &JxW = feFaceValues.get_JxW_values();
	const vector<dealii::Tensor<1, dim> > &normals =
			feFaceValues.get_all_normal_vectors();

	dealii::FullMatrix<double> anotherDirectionCellFaceMatrix(
			feFaceValues.dofs_per_cell);

	for (size_t beta = 1; beta < stencil.getQ(); beta++) {

		anotherDirectionCellFaceMatrix = 0;

		// calculate prefactor
		for (size_t q = 0; q < feFaceValues.n_quadrature_points; q++) {
			size_t thisDoF = q_index_to_facedof.at(q);

			// get velocity at boundary point
			dealii::Vector<double> velocity(dim);
			DirichletBoundary<dim>::getBoundaryVelocity()->vector_value(
					feFaceValues.quadrature_point(q), velocity);

			// calculate matrix entries
			double prefactor = JxW.at(q);
			double exn = 0.0;
			// calculate scalar product
			for (size_t i = 0; i < dim; i++) {	// TODO efficient multiplication
				exn += normals.at(q)[i] * stencil.getDirection(alpha)(i);
			}
			prefactor *= exn;

			double ea_x_u = stencil.getDirection(alpha)
					* velocity;

			if (useCentralFlux) {
				assert(false);
			} else if (exn < 0) {
				anotherDirectionCellFaceMatrix(thisDoF, thisDoF) = -prefactor
						* 2 * stencil.getWeight(alpha) * ea_x_u
						/ stencil.getSpeedOfSoundSquare();
				if (beta == stencil.getIndexOfOppositeDirection(alpha))
					anotherDirectionCellFaceMatrix(thisDoF, thisDoF) = anotherDirectionCellFaceMatrix(thisDoF, thisDoF) - prefactor;
				if (beta == alpha)
					anotherDirectionCellFaceMatrix(thisDoF, thisDoF) = anotherDirectionCellFaceMatrix(thisDoF, thisDoF) + prefactor;
			}
		}

		vector<dealii::types::global_dof_index> localDoFIndices(
				feFaceValues.dofs_per_cell);
		cell->get_dof_indices(localDoFIndices);

		/// Distribute to global matrix
		for (size_t i = 0; i < feFaceValues.dofs_per_cell; i++) {
			if (cell->get_fe().has_support_on_face(i, faceNumber)) {
				systemMatrix.block(alpha - 1,
						beta - 1).add(
						localDoFIndices[i], localDoFIndices[i],
						anotherDirectionCellFaceMatrix(i, i)
								* inverseLocalMassMatrix.at(i));
			}
		}
	} /* for beta */
}


// Explicit instantiation
template class DirichletBoundaryU<2> ;
template class DirichletBoundaryU<3> ;

} /* namespace natrium */
