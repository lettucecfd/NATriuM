/*
 * DirichletBoundaryRho.cpp
 *
 *  Created on: 08.12.2015
 *      Author: akraem3m
 */

#include "DirichletBoundaryRho.h"

namespace natrium {

template<size_t dim>
DirichletBoundaryRho<dim>::DirichletBoundaryRho(size_t boundaryIndicator,
		boost::shared_ptr<dealii::Function<dim> > boundaryDensity) :
		DirichletBoundary<dim>(boundaryIndicator, boundaryDensity,
		NULL, BoundaryTools::COUPLE_ALL_DISTRIBUTIONS,
				BoundaryTools::COUPLE_ONLY_SINGLE_POINTS) {

}

/// constructor
template<size_t dim>
DirichletBoundaryRho<dim>::DirichletBoundaryRho(size_t boundaryIndicator,
		double rho) :
		DirichletBoundaryRho<dim>(boundaryIndicator,
				boost::make_shared<BoundaryDensity<dim> >(rho)) {

}

template<size_t dim>
DirichletBoundaryRho<dim>::~DirichletBoundaryRho() {
}

template<size_t dim> void DirichletBoundaryRho<dim>::assembleBoundary(
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

	dealii::FullMatrix<double> cellFaceMatrix(feFaceValues.dofs_per_cell);
	dealii::FullMatrix<double> anotherDirectionCellFaceMatrix(
			feFaceValues.dofs_per_cell);

	for (size_t beta = 0; beta < stencil.getQ(); beta++) {

		double ea_x_eb = stencil.getDirection(alpha)
				* stencil.getDirection(beta);
		if (fabs(ea_x_eb) < 1e-10)
			continue;

		cellFaceMatrix = 0;
		anotherDirectionCellFaceMatrix = 0;

		// calculate prefactor
		for (size_t q = 0; q < feFaceValues.n_quadrature_points; q++) {
			size_t thisDoF = q_index_to_facedof.at(q);

			// get density at boundary point
			double density =
					DirichletBoundary<dim>::getBoundaryDensity()->value(
							feFaceValues.quadrature_point(q));

			// calculate matrix entries
			double prefactor = JxW.at(q);
			double exn = 0.0;
			// calculate scalar product
			for (size_t i = 0; i < dim; i++) {	// TODO efficient multiplication
				exn += normals.at(q)[i] * stencil.getDirection(alpha)(i);
			}
			prefactor *= exn;

			if (useCentralFlux) {
				// not implemented (also not needed)
				assert(false);
			} else if (exn < 0) {
				cellFaceMatrix(thisDoF, thisDoF) = prefactor * 2
						* stencil.getWeight(alpha) * density * ea_x_eb
						/ stencil.getSpeedOfSoundSquare();
				anotherDirectionCellFaceMatrix(thisDoF, thisDoF) = -prefactor
						* 2 * stencil.getWeight(alpha) * density * ea_x_eb
						/ stencil.getSpeedOfSoundSquare();
			}
		}

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
						anotherDirectionCellFaceMatrix(i, i)
								* inverseLocalMassMatrix.at(i));
			}
		}
	} /* for beta */
}

// Explicit instantiation
template class DirichletBoundaryRho<2> ;
template class DirichletBoundaryRho<3> ;

} /* namespace natrium */
