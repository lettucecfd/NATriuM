/*
 * LinearBoundaryRhoU.cpp
 *
 *  Created on: 08.12.2015
 *      Author: akraem3m
 */

#include "VelocityNeqBounceBack.h"
#include "BoundaryFlags.h"
#include "BoundaryTools.h"

namespace natrium {

template<size_t dim>
VelocityNeqBounceBack<dim>::VelocityNeqBounceBack(size_t boundaryIndicator,
		boost::shared_ptr<dealii::Function<dim> > boundaryVelocity) :
		Boundary<dim>(boundaryIndicator, NONCONSTANT_VELOCITY_NEQ_BB,
				PrescribedBoundaryValues<dim>(boundaryVelocity) ) {

	//assert(not Boundary<dim>::getBoundaryValues().getPressure());
	assert(Boundary<dim>::getBoundaryValues().getVelocity());

}

/// constructor
template<size_t dim>
VelocityNeqBounceBack<dim>::VelocityNeqBounceBack(size_t boundaryIndicator,
		const dealii::Vector<double>& velocity) :
		VelocityNeqBounceBack<dim>(boundaryIndicator,
				boost::make_shared<BoundaryTools::BoundaryVelocity<dim> >(velocity)) {

}

template<size_t dim>
VelocityNeqBounceBack<dim>::VelocityNeqBounceBack(size_t boundaryIndicator,
		const dealii::Tensor<1,dim>& velocity):
	VelocityNeqBounceBack<dim>(boundaryIndicator,
					boost::make_shared<BoundaryTools::BoundaryVelocity<dim> >(velocity))  {
}


template<size_t dim>
VelocityNeqBounceBack<dim>::~VelocityNeqBounceBack() {
}

template<size_t dim> void VelocityNeqBounceBack<dim>::assembleBoundary(
		size_t alpha,
		const typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
		size_t faceNumber, dealii::FEFaceValues<dim>& feFaceValues,
		const Stencil& stencil,
		const std::map<size_t, size_t>& q_index_to_facedof,
		const vector<double> & inverseLocalMassMatrix,
		distributed_sparse_block_matrix& systemMatrix,
		distributed_block_vector& systemVector, bool useCentralFlux) {
	// let the feFaceValues object calculate all the values needed at the boundary
	feFaceValues.reinit(cell, faceNumber);
	const vector<double> &JxW = feFaceValues.get_JxW_values();
	const vector<dealii::Tensor<1, dim> > &normals =
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
		double density = 1.0; // LinearFluxBoundary<dim>::getBoundaryDensity()->value(
				// feFaceValues.quadrature_point(q));
		dealii::Vector<double> velocity(dim);

		assert(Boundary<dim>::getBoundaryValues().getVelocity());
		Boundary<dim>::getBoundaryValues().getVelocity()->vector_value(
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


template <size_t dim>
void VelocityNeqBounceBack<dim>::calculateBoundaryValues(
		FEBoundaryValues<dim>& fe_boundary_values, size_t q_point,
		const LagrangianPathDestination& destination, double eps,
		double t) {

	const GlobalBoundaryData& data = fe_boundary_values.getData();
	const Stencil& stencil = data.m_stencil;

	// density unity
	double rho = 1;

	// get velocity
	dealii::Tensor<1, dim> u;
	// set time of this function object for time-dependent boundary conditions


	assert(Boundary<dim>::getBoundaryValues().getVelocity());
	Boundary<dim>::getBoundaryValues().getVelocity()->set_time(
			t - eps);
	// evaluate boundary conditions
	dealii::Vector<double> tmp(dim);
	Boundary<dim>::getBoundaryValues().getVelocity()->vector_value(
			fe_boundary_values.getPoint(q_point), tmp);


	// vector to tensor
	for (size_t i = 0; i < dim; i++) {
		u[i] = tmp(i);
	}


	// calculate vector entry
	double exu = 0.0;
	double uxu = 0.0;
	double cs2 = stencil.getSpeedOfSoundSquare();

	// calculate scalar product
	for (size_t i = 0; i < dim; i++) {	// TODO efficient multiplication
		exu += u[i] * stencil.getDirection(destination.direction)(i);
		uxu += u[i] * u[i];
	}


	double to_add = 2 * stencil.getWeight(destination.direction)
			* rho * exu / stencil.getSpeedOfSoundSquare();


	 //equilibrium for testing:
	  fe_boundary_values.getData().m_fnew.at(destination.direction)(
			destination.index) = stencil.getWeight(destination.direction) * rho *
			(1 + exu / cs2 + u*u / (2 * cs2) +  (exu * exu) / (2*cs2*cs2));
		/*
    if (not destination.domain_corner) {
        fe_boundary_values.getData().m_fnew.at(destination.direction)(
                destination.index) =
                fe_boundary_values.getData().m_fnew.at(destination.direction)(
                        destination.index) + to_add;
    } */

}

// Explicit instantiation
template class VelocityNeqBounceBack<2> ;
template class VelocityNeqBounceBack<3> ;

} /* namespace natrium */
