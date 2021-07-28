/*
 * BoundaryFlags.h
 *
 *  Created on: 29.11.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_BOUNDARIES_BOUNDARYFLAGS_H_
#define LIBRARY_NATRIUM_BOUNDARIES_BOUNDARYFLAGS_H_

#include "BoundaryTools.h"

namespace natrium {

/**
 * @short An enum that is used to flag values at the boundary.
 * 			BoundaryFlags is used by semi-Lagrangian boundary conditions to flag prescribed values
 * 			or values that are to be calculated from the flow field.
 */
enum BoundaryFlags {
	only_distributions = 0x0000,
	boundary_rho = 0x0001,
	boundary_u = 0x0002,
	boundary_du_dx = 0x0004,
	boundary_drho_dt = 0x0008,
	boundary_du_dt = 0x0010,
	boundary_d2u_dxdt = 0x0020,
	boundary_p = 0x0040,
	boundary_dp_dt = 0x0080,
	boundary_pUT = 0x0100
};

inline BoundaryFlags
operator | (BoundaryFlags f1, BoundaryFlags f2){
	return static_cast<BoundaryFlags> (
			static_cast<unsigned int> (f1)  |
			static_cast<unsigned int> (f2));
}
inline BoundaryFlags&
operator |= (BoundaryFlags& f1, BoundaryFlags f2){
	f1 = f1 | f2;
	return f1;
}
inline BoundaryFlags
operator & (BoundaryFlags f1, BoundaryFlags f2){
	return static_cast<BoundaryFlags> (
			static_cast<unsigned int> (f1)  &
			static_cast<unsigned int> (f2));
}
inline BoundaryFlags &
operator &= (BoundaryFlags &f1, BoundaryFlags f2){
	f1 = f1 & f2;
	return f1;
}

/**
 * @short A class that contains the values that are prescribed by a boundary condition
 */
template<size_t dim>
class PrescribedBoundaryValues {
private:
	BoundaryFlags m_prescribedValues;
	boost::shared_ptr<dealii::Function<dim> > m_pressure;
	boost::shared_ptr<dealii::Function<dim> > m_velocity;
	boost::shared_ptr<dealii::Function<dim> > m_temperature;
	boost::shared_ptr<dealii::TensorFunction<2, dim> > m_velocityGradient;

public:
	/**
	 * @short Empty constructor. No values to be obtained (used e.g. for first-order bounce back)
	 */
	PrescribedBoundaryValues() {
		m_prescribedValues = only_distributions;
	}

	/**
	 * @short Copy constructor
	 */
	PrescribedBoundaryValues(const PrescribedBoundaryValues<dim>& other) {
		m_prescribedValues = other.m_prescribedValues;
		m_pressure = other.m_pressure;
		m_velocity = other.m_velocity;
		m_velocityGradient = other.m_velocityGradient;
	}

	/**
	 * @short Constructor using a function.
	 * If the number of components in the function is 1, it is interpreted as a pressure boundary condition.
	 * If the number is dim, it is interpreted as a velocity boundary condition.
	 */
	PrescribedBoundaryValues(boost::shared_ptr<dealii::Function<dim> > function) {
		m_prescribedValues = only_distributions;
		if (function->n_components == 1) {
			m_pressure = function;
			m_prescribedValues = boundary_p;
		} else if (function->n_components == dim) {
			m_velocity = function;
			m_prescribedValues = boundary_u;
		} else {
			throw NATriuMException(
					"The function that was to prescribe quantities at the boundary had the wrong "
							"number of components (has to be either 1 or dim).");
		}
	}

	/**
	 * @short Constructor using a tensor function defining the Jacobian du/dx at the boundary.
	 */
	PrescribedBoundaryValues(
			boost::shared_ptr<dealii::TensorFunction<2, dim> > function) {
		m_prescribedValues = boundary_du_dx;
		m_velocityGradient = function;
	}

	/**
	 * @short Constructor using a double, i.e. the pressure
	 * @note The function instance is created by using BoundaryTools::BoundaryPressure
	 */
	PrescribedBoundaryValues(double p) {
		m_prescribedValues = boundary_p;
		m_pressure = boost::make_shared<BoundaryTools::BoundaryPressure<dim> >(
				p);
	}

	/**
	 * @short Constructor using a tensor, i.e. the velocity.
	 * @note The function instance is created by using BoundaryTools::BoundaryVelocity
	 */
	PrescribedBoundaryValues(const dealii::Tensor<1, dim>& u) {
		m_prescribedValues = boundary_u;
		m_velocity = boost::make_shared<BoundaryTools::BoundaryVelocity<dim> >(
				u);
	}

    /**
 * @short Constructor using a tensor, i.e. the velocity.
 * @note The function instance is created by using BoundaryTools::BoundaryVelocity
 */
    PrescribedBoundaryValues(double p, const dealii::Tensor<1, dim>& u, double T) {
        m_prescribedValues = boundary_pUT;
        m_pressure = boost::make_shared<BoundaryTools::BoundaryPressure<dim> >(
                p);
        m_velocity = boost::make_shared<BoundaryTools::BoundaryVelocity<dim> >(
                u);
        m_temperature = boost::make_shared<BoundaryTools::BoundaryTemperature<dim> >(
                p);
    }

	/**
	 * @short Constructor using a dealii::Vector, i.e. the velocity.
	 * @note The function instance is created by using BoundaryTools::BoundaryVelocity
	 */
	PrescribedBoundaryValues(const dealii::Vector<double>& u) {
		m_prescribedValues = boundary_u;
		m_velocity = boost::make_shared<BoundaryTools::BoundaryVelocity<dim> >(
				u);
	}

	/**
	 * @short Constructor using a tensor, i.e. the Jacobian
	 * @note The function instance is created by using BoundaryTools::BoundaryVelocityGradient
	 */
	PrescribedBoundaryValues(const dealii::Tensor<2, dim>& dudx) {
		m_prescribedValues = boundary_du_dx;
		m_velocityGradient = boost::make_shared<
				BoundaryTools::BoundaryVelocityGradient<dim> >(dudx);
	}

	/**
	 * @short Adding prescribed values at corner nodes. Not implemented yet.
	 */
	PrescribedBoundaryValues& operator+(const PrescribedBoundaryValues& other) {
		PrescribedBoundaryValues result;
		result.m_prescribedValues = this->m_prescribedValues
				| other.m_prescribedValues;
		LOG(ERROR) << "Prescribed Values operator+ is not yet complete." << endl;
		// TODO combine values
		return result;
	}

	/**
	 * @short get the boundary flags that define the prescribed values
	 */
	BoundaryFlags getPrescribedValues() const {
		return m_prescribedValues;
	}

	/**
	 * @short get the function that defines the prescribed pressure
	 */
	boost::shared_ptr<dealii::Function<dim> >& getPressure()  {
		return m_pressure;
	}

	/**
	 * @short get the function that defines the prescribed velocity
	 */
	boost::shared_ptr<dealii::Function<dim> >& getVelocity()  {
		return m_velocity;
	}
/**
	 * @short get the function that defines the prescribed temperature
	 */
    boost::shared_ptr<dealii::Function<dim> >& getTemperature()  {
        return m_temperature;
    }

	/**
	 * @short get the function that defines the prescribed Jacobian du/dx
	 */
	boost::shared_ptr<dealii::TensorFunction<2, dim> >& getVelocityGradient()  {
		return m_velocityGradient;
	}
};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_BOUNDARIES_BOUNDARYFLAGS_H_ */
