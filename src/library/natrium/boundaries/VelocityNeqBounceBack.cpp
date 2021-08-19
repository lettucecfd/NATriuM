/*
 * LinearBoundaryRhoU.cpp
 *
 *  Created on: 08.12.2015
 *      Author: akraem3m
 */

#include "VelocityNeqBounceBack.h"
#include "BoundaryFlags.h"
#include "BoundaryTools.h"
#include "../collision_advanced/AuxiliaryCollisionFunctions.h"

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

	const double temperature = 1.0;
	const double density = rho;
	double feq;
    const double weight =  stencil.getWeight(destination.direction);
    std::array<double,dim> e_vel ={0.0};
    for (size_t i = 0; i < dim; i++) {
       e_vel[i] = stencil.getDirection(destination.direction)(i);
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


    const std::array<std::array<size_t,dim>, dim> eye = unity_matrix<dim>();
    double uu_term = 0.0;
    for (size_t j = 0; j < dim; j++) {
        uu_term += -(u[j] * u[j])
                   / (2.0 * cs2);
    }



        double T1 = cs2*(temperature-1);

        double ue_term = 0.0;
        for (size_t j = 0; j < dim; j++) {
            ue_term += (u[j] * e_vel[j]) / cs2;
        }
        feq = weight * density * (1 + ue_term * (1 + 0.5 * (ue_term)) + uu_term);
        for (size_t alp = 0; alp < dim; alp++){
            for (size_t bet = 0; bet < dim; bet++){
                feq+=density*weight/(2.0*cs2)*((temperature-1)*eye[alp][bet]*e_vel[alp]*e_vel[bet]-cs2*eye[alp][bet]*(temperature-1));
                for (size_t gam = 0; gam < dim; gam++){

                    feq += weight * density / (6. * cs2 * cs2 * cs2) *
                              (u[alp] * u[bet] * u[gam]
                               + T1 *
                                 (eye[alp][bet] * u[gam] + eye[bet][gam] * u[alp] +
                                  eye[alp][gam] * u[bet])) * (e_vel[alp] * e_vel[bet] * e_vel[gam] - cs2 *
                                                                                                         (e_vel[gam] *
                                                                                                          eye[alp][bet] +
                                                                                                          e_vel[bet] *
                                                                                                          eye[alp][gam] +
                                                                                                          e_vel[alp] *
                                                                                                          eye[bet][gam]));

                    for (size_t det = 0; det < dim; det++)
                    {
                        double power4 = e_vel[alp]*e_vel[bet]*e_vel[gam]*e_vel[det];
                        double power2 = e_vel[alp]*e_vel[bet]*eye[gam][det]
                                        +e_vel[alp]*e_vel[gam]*eye[bet][det]
                                        +e_vel[alp]*e_vel[det]*eye[bet][gam]
                                        +e_vel[bet]*e_vel[gam]*eye[alp][det]
                                        +e_vel[bet]*e_vel[det]*eye[alp][gam]
                                        +e_vel[gam]*e_vel[det]*eye[alp][bet];
                        double power0 = eye[alp][bet]*eye[gam][det]+eye[alp][gam]*eye[bet][det]+eye[alp][det]*eye[bet][gam];
                        double u4    = u[alp]*u[bet]*u[gam]*u[det];
                        double u2 = u[alp]*u[bet]*eye[gam][det]+u[alp]*u[gam]*eye[bet][det]+u[alp]*u[det]*eye[bet][gam]+u[bet]*u[gam]*eye[alp][det]+u[bet]*u[det]*eye[alp][gam]+u[gam]*u[det]*eye[alp][bet];
                        double multieye= eye[alp][bet]*eye[gam][det]+eye[alp][gam]*eye[bet][det]+eye[alp][det]*eye[bet][gam];

                        feq+= weight * density /(24.*cs2*cs2*cs2*cs2)*(power4-cs2*power2+cs2*cs2*power0)*(u4+T1*(u2+T1*multieye));


                    }

                }
            }

        }


	 //equilibrium for testing:
	 if (uxu>0.000001)
	  fe_boundary_values.getData().m_fnew.at(destination.direction)(
			destination.index) = feq; //stencil.getWeight(destination.direction) * rho *
			//(1 + exu / cs2 + u*u / (2 * cs2) +  (exu * exu) / (2*cs2*cs2));
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
