/*
 * SLEquilibriumBoundary.h
 *
 *  Created on: 28.11.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_BOUNDARIES_SLEQUILIBRIUMBOUNDARY_H_
#define LIBRARY_NATRIUM_BOUNDARIES_SLEQUILIBRIUMBOUNDARY_H_

#include "../utilities/BasicNames.h"
#include "BoundaryFlags.h"
#include "../collision_advanced/AuxiliaryCollisionFunctions.h"


namespace natrium {

/**
 * @short Boundary condition that sets an equilibrium with respect to a prescribed velocity or pressure.
 */
template<size_t dim>
class SLEquilibriumBoundary: public Boundary<dim> {

public:
	SLEquilibriumBoundary(size_t boundary_id, dealii::Tensor<1, dim>& velocity, double temperature = 1.0) :
			Boundary<dim>(boundary_id, VELOCITY_EQUILIBRIUM_BOUNDARY, PrescribedBoundaryValues<dim>(velocity)),
			        m_temperature(temperature) {
	}

    SLEquilibriumBoundary(size_t boundary_id, const dealii::Vector<double>& velocity, double temperature = 1.0) :
            Boundary<dim>(boundary_id, VELOCITY_EQUILIBRIUM_BOUNDARY, PrescribedBoundaryValues<dim>(velocity)),
                    m_temperature(temperature) {
    }

	//virtual ~SLEquilibriumBoundary();

	/////////////////////////////////////////////////
	// FLAGS ////////////////////////////////////////
	/////////////////////////////////////////////////

	/** @short is the boundary a periodic boundary ?
	 */
	virtual bool isPeriodic() const {
		return false;
	}

	/** @short is the boundary a linear flux boundary as in SEDG-LBM (affine linear in the distributions f)?
	 */
	virtual bool isDGSupported() const {
		return false;
	}

	/** @short is the boundary set up to work with semi-Lagrangian streaming
	 */
	virtual bool isSLSupported() const {
		return true;
	}



    virtual BoundaryFlags getUpdateFlags() const {
        BoundaryFlags flags = only_distributions;
        return flags;
    }

	virtual void calculateBoundaryValues(FEBoundaryValues<dim>& fe_boundary_values, size_t q_point,
                                         const LagrangianPathDestination& destination, double eps, double t) {

        const GlobalBoundaryData& data = fe_boundary_values.getData();
        const Stencil& stencil = data.m_stencil;

        const double scaling = stencil.getScaling();

        // get velocity
        dealii::Tensor<1, dim> u;
        assert(Boundary<dim>::getBoundaryValues().getVelocity());
        // set time of this function object for time-dependent boundary conditions
        Boundary<dim>::getBoundaryValues().getVelocity()->set_time(t - eps);
        // evaluate boundary conditions
        dealii::Vector<double> tmp(dim);
        Boundary<dim>::getBoundaryValues().getVelocity()->vector_value(fe_boundary_values.getPoint(q_point), tmp);
        // vector to tensor
        for (size_t i = 0; i < dim; i++) {
            u[i] = tmp(i);
        }
        // density unity
        double rho = 1;
        const double temperature = m_temperature;
        const double density = rho;
        double feq;
        const double weight = stencil.getWeight(destination.direction);
        std::array<double, dim> e_vel = {};
        for (size_t i = 0; i < dim; i++) {
            e_vel[i] = stencil.getDirection(destination.direction)(i);
        }
        double cs2 = stencil.getSpeedOfSoundSquare();

            const array<array<size_t, dim>, dim> eye = unity_matrix<dim>();
            double uu_term = 0.0;
            for (size_t j = 0; j < dim; j++) {
                uu_term += -(u[j] * u[j]) / (2.0 * cs2);
            }

            double T1 = cs2 * (temperature - 1);

            double ue_term = 0.0;
            for (size_t j = 0; j < dim; j++) {
                ue_term += (u[j] * e_vel[j]) / cs2;
            }
            feq = weight * density * (1 + ue_term * (1 + 0.5 * (ue_term)) + uu_term);
            for (size_t alp = 0; alp < dim; alp++) {
                for (size_t bet = 0; bet < dim; bet++) {
                    feq += density * weight / (2.0 * cs2) * (temperature - 1)
                            * (eye[alp][bet] * e_vel[alp] * e_vel[bet] - cs2 * eye[alp][bet]);
                    for (size_t gam = 0; gam < dim; gam++) {
                        feq += weight * density / (6. * cs2 * cs2 * cs2) *
                                (u[alp] * u[bet] * u[gam]
                                    + T1 * (eye[alp][bet] * u[gam]
                                          + eye[bet][gam] * u[alp]
                                          + eye[alp][gam] * u[bet])) *
                                (e_vel[alp] * e_vel[bet] * e_vel[gam]
                                    - cs2 * (e_vel[gam] * eye[alp][bet]
                                           + e_vel[bet] * eye[alp][gam]
                                           + e_vel[alp] * eye[bet][gam]));

                        for (size_t det = 0; det < dim; det++) {
                            double power4 = e_vel[alp] * e_vel[bet] * e_vel[gam] * e_vel[det];
                            double power2 = e_vel[alp] * e_vel[bet] * eye[gam][det]
                                            + e_vel[alp] * e_vel[gam] * eye[bet][det]
                                            + e_vel[alp] * e_vel[det] * eye[bet][gam]
                                            + e_vel[bet] * e_vel[gam] * eye[alp][det]
                                            + e_vel[bet] * e_vel[det] * eye[alp][gam]
                                            + e_vel[gam] * e_vel[det] * eye[alp][bet];
                            double power0 = eye[alp][bet] * eye[gam][det]
                                          + eye[alp][gam] * eye[bet][det]
                                          + eye[alp][det] * eye[bet][gam];
                            double u4 = u[alp] * u[bet] * u[gam] * u[det];
                            double u2 = u[alp] * u[bet] * eye[gam][det] + u[alp] * u[gam] * eye[bet][det] +
                                        u[alp] * u[det] * eye[bet][gam] + u[bet] * u[gam] * eye[alp][det] +
                                        u[bet] * u[det] * eye[alp][gam] + u[gam] * u[det] * eye[alp][bet];
                            double multieye = eye[alp][bet] * eye[gam][det]
                                            + eye[alp][gam] * eye[bet][det]
                                            + eye[alp][det] * eye[bet][gam];

                            feq += weight * density / (24. * cs2 * cs2 * cs2 * cs2) *
                                   (power4 - cs2 * power2 + cs2 * cs2 * power0) * (u4 + T1 * (u2 + T1 * multieye));
                        }
                    }
                }
            }
            fe_boundary_values.getData().m_fnew.at(destination.direction)(destination.index) = feq;
            //// NEW: calculate geq
            const double gamma = 1.4;
            const double C_v = 1. / (gamma - 1.0);
            fe_boundary_values.getData().m_g.at(destination.direction)(destination.index) =
                    feq * temperature * (2.0 * C_v - dim);
    }
private:
    double m_temperature;
};

} /* namespace natrium */


#endif /* LIBRARY_NATRIUM_BOUNDARIES_SLEQUILIBRIUMBOUNDARY_H_ */
