/*
 * SLBoundary.h
 *
 *  Created on: 19.09.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_BOUNDARIES_SLBOUNDARY_H_
#define LIBRARY_NATRIUM_BOUNDARIES_SLBOUNDARY_H_

#include "Boundary.h"

#include "deal.II/base/function.h"
#include "deal.II/base/tensor_function.h"

#include "../solver/DistributionFunctions.h"
#include "../stencils/Stencil.h"
#include "../advection/SemiLagrangianTools.h"
#include "../utilities/BasicNames.h"
#include "../advection/AdvectionOperator.h"

namespace natrium {

template<size_t dim>
struct PrescribedQuantities {
public:
	bool prescribedPressure;
	bool prescribedVelocity;
	bool prescribedVelocityGradient;
	boost::shared_ptr<dealii::Function<dim> > pressure;
	boost::shared_ptr<dealii::Function<dim> > velocity;
	boost::shared_ptr<dealii::TensorFunction<2, dim> > velocityGradient;

	PrescribedQuantities(const PrescribedQuantities<dim>& other) {
		prescribedPressure = other.prescribedPressure;
		prescribedVelocity = other.prescribedVelocity;
		prescribedVelocityGradient = other.prescribedVelocityGradient;
		pressure = other.pressure;
		velocity = other.velocity;
		velocityGradient = other.velocityGradient;
	}
	PrescribedQuantities(boost::shared_ptr<dealii::Function<dim> > function) {
		prescribedPressure = false;
		prescribedVelocity = false;
		prescribedVelocityGradient = false;
		if (function->n_components == 1) {
			pressure = function;
			prescribedPressure = true;
		} else if (function->n_components == dim) {
			velocity = function;
			prescribedVelocity = true;
		} else {
			throw NATriuMException(
					"The function that was to prescribe quantities at the boundary had the wrong "
							"number of components (has to be either 1 or dim).");
		}
	}

	PrescribedQuantities(
			boost::shared_ptr<dealii::TensorFunction<2, dim> > function) {
		prescribedPressure = false;
		prescribedVelocity = false;
		prescribedVelocityGradient = true;
		velocityGradient = function;
	}

	PrescribedQuantities(double p) {
		prescribedPressure = true;
		prescribedVelocity = false;
		prescribedVelocityGradient = false;
		pressure = boost::make_shared<BoundaryTools::BoundaryPressure<dim> >(p);
	}

	PrescribedQuantities(const dealii::Tensor<1, dim>& u) {
		prescribedPressure = false;
		prescribedVelocity = true;
		prescribedVelocityGradient = false;
		velocity = boost::make_shared<BoundaryTools::BoundaryVelocity<dim> >(u);
	}

	PrescribedQuantities(const dealii::Vector<double>& u) {
		assert(u.size() == dim);
		prescribedPressure = false;
		prescribedVelocity = true;
		prescribedVelocityGradient = false;
		velocity = boost::make_shared<BoundaryTools::BoundaryVelocity<dim> >(u);
	}

	PrescribedQuantities(const dealii::Tensor<2, dim>& dudx) {
		prescribedPressure = false;
		prescribedVelocity = false;
		prescribedVelocityGradient = true;
		velocityGradient = boost::make_shared<
				BoundaryTools::BoundaryVelocityGradient<dim> >(dudx);
	}
};

/**
 * @short  variables that are required in each iteration
 */
template<size_t dim>
struct LocalBoundaryData {
	// constants
	size_t dofs_per_cell;
	std::vector<dealii::types::global_dof_index>  local_dofs;

	// macroscopic variables and dof variables
	double rho;
	dealii::Tensor<1, dim> rho_u;
	dealii::Tensor<1, dim> u;
	vector<double> rho_dof;
	vector<dealii::Tensor<1, dim> > rho_u_dof;
	vector<dealii::Tensor<1, dim> > u_dof;

	// wall values and temporal derivatives
	double rho_w;
	double drhodt_w;
	dealii::Tensor<1, dim> u_w;
	dealii::Tensor<1, dim> dudt_w;
	dealii::Tensor<2, dim> dudx_w;
	dealii::Tensor<2, dim> d2udxdt_w;

	LocalBoundaryData() {
		dofs_per_cell = 0;

		// macroscopic variables and dof variables
		rho = 0;
		rho_u = 0;
		u = 0;
		/*rho_dof.clear();
		rho_u_dof.clear();
		u_dof.clear();
		local_dofs.clear();
		 */

		// wall values and temporal derivatives
		rho_w = 0;
		drhodt_w = 0;
		u_w = 0;
		dudt_w = 0;
		dudx_w = 0;
		d2udxdt_w = 0;
	}

	virtual ~LocalBoundaryData() {
	}

	void resize(size_t n) {
		rho_dof.resize(n);
		rho_u_dof.resize(n);
		u_dof.resize(n);
		local_dofs.resize(n);

	}
	void setToZero() {
		// macroscopic variables and dof variables
		rho = 0;
		rho_u = 0;
		u = 0;

		// wall values and temporal derivatives
		rho_w = 0;
		drhodt_w = 0;
		u_w = 0;
		dudt_w = 0;
		dudx_w = 0;
		d2udxdt_w = 0;
	}

	void clear() {
		dofs_per_cell = 0;

		// macroscopic variables and dof variables
		rho = 0;
		rho_u = 0;
		u = 0;
		rho_dof.clear();
		rho_u_dof.clear();
		u_dof.clear();
		local_dofs.clear();

		// wall values and temporal derivatives
		rho_w = 0;
		drhodt_w = 0;
		u_w = 0;
		dudt_w = 0;
		dudx_w = 0;
		d2udxdt_w = 0;
	}
};

struct GlobalBoundaryData {
	const DistributionFunctions& m_fold;
	DistributionFunctions& m_fnew;
	const Stencil& m_stencil;
	double m_viscosity;
	double m_dt;
	size_t m_Q;
	double m_cs2;
	GlobalBoundaryData(const DistributionFunctions& f_old,
			DistributionFunctions& f_new, const Stencil& stencil,
			double viscosity, double dt) :
			m_fold(f_old), m_fnew(f_new), m_stencil(stencil) {
		m_viscosity = viscosity;
		m_dt = dt;
		m_cs2 = stencil.getSpeedOfSoundSquare();
		m_Q = stencil.getQ();
	}
	virtual ~GlobalBoundaryData() {

	}
};

template<size_t dim>
class SLBoundary: public Boundary<dim> {
private:
	size_t m_boundaryIndicator;
	PrescribedQuantities<dim> m_prescribedQuantities;
public:

	SLBoundary(size_t boundaryIndicator,
			const PrescribedQuantities<dim>& quantities) :
			m_boundaryIndicator(boundaryIndicator), m_prescribedQuantities(
					quantities) {

	}

	virtual ~SLBoundary() {
	}
	;

	virtual bool isPeriodic() const {
		return false;
	}
	virtual bool isLinearFluxBoundary() const {
		return false;
	}
	virtual bool isSLBoundary() const {
		return true;
	}

	size_t getBoundaryIndicator() const {
		return m_boundaryIndicator;
	}

	virtual void calculateBoundaryValues(const GlobalBoundaryData& g,
				LocalBoundaryData<dim>& b, const dealii::FEValues<dim>& fe_values,
				size_t q_point, const LagrangianPathDestination& destination,
				double dt) = 0;

	virtual dealii::UpdateFlags getUpdateFlags() const = 0;

	const PrescribedQuantities<dim>& getPrescribedQuantities() const {
		return m_prescribedQuantities;
	}
};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_BOUNDARIES_SLBOUNDARY_H_ */

