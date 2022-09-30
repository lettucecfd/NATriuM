/*
 * ThermalBounceBack.cpp
 *
 *  Created on: 29.10.2022
 *      Author: Dominik Wilde
 */

#include "ThermalBounceBack.h"
#include "BoundaryFlags.h"
#include "BoundaryTools.h"
#include "../collision_advanced/AuxiliaryCollisionFunctions.h"
#include "../collision_advanced/Equilibria.h"

namespace natrium {

template<size_t dim>
ThermalBounceBack<dim>::ThermalBounceBack(size_t boundaryIndicator,
		boost::shared_ptr<dealii::Function<dim> > boundaryVelocity) :
		Boundary<dim>(boundaryIndicator, THERMAL_BB,
				PrescribedBoundaryValues<dim>(boundaryVelocity) ) {

	//assert(not Boundary<dim>::getBoundaryValues().getPressure());
	assert(Boundary<dim>::getBoundaryValues().getVelocity());

}

/// constructor
template<size_t dim>
ThermalBounceBack<dim>::ThermalBounceBack(size_t boundaryIndicator,
		const dealii::Vector<double>& velocity) :
		ThermalBounceBack<dim>(boundaryIndicator,
				boost::make_shared<BoundaryTools::BoundaryVelocity<dim> >(velocity)) {

}

template<size_t dim>
ThermalBounceBack<dim>::ThermalBounceBack(size_t boundaryIndicator,
		const dealii::Tensor<1,dim>& velocity):
	ThermalBounceBack<dim>(boundaryIndicator,
					boost::make_shared<BoundaryTools::BoundaryVelocity<dim> >(velocity))  {
}


template<size_t dim>
ThermalBounceBack<dim>::~ThermalBounceBack() {
}


template <size_t dim>
void ThermalBounceBack<dim>::calculateBoundaryValues(
		FEBoundaryValues<dim>& fe_boundary_values, size_t q_point,
		const LagrangianPathDestination& destination, double eps,
		double t) {

	const GlobalBoundaryData& data = fe_boundary_values.getData();
	const Stencil& stencil = data.m_stencil;
    const double scaling = stencil.getScaling();
    const double cs2 = stencil.getSpeedOfSoundSquare() / (scaling * scaling);
    const double gamma = 1.4;
    assert(stencil.getQ()==45);
    std::array<double,45> f_destination, g_destination, feq,w ;
    for (int i=0; i<45; i++) {
        f_destination[i] = fe_boundary_values.getData().m_fnew.at(i)(destination.index);
        g_destination[i] = fe_boundary_values.getData().m_g.at(i)(destination.index);
        w[i]=stencil.getWeight(i);
    }

    const double rho = calculateDensity<45>(f_destination);
    std::array<double,dim> u_local;
    std::array<std::array<double,dim>,45> e = getParticleVelocitiesWithoutScaling<dim,45>(stencil);
    calculateVelocity<dim,45>(f_destination,u_local,rho,e);

    const double T_local = calculateTemperature<dim,45>(f_destination,g_destination,u_local,rho,e,cs2,1.4);

    QuarticEquilibrium<dim,45> eq(cs2,e);
    eq.polynomial(feq,rho,u_local,T_local,e,w,cs2);
    for (int i=0; i<45; i++) {
        f_destination[i] -= feq[i];
    }
    eq.polynomial(feq,rho,u_local,1.0,e,w,cs2);
    for (int i=0; i<45; i++) {
        //f_destination[i] += feq[i];
        fe_boundary_values.getData().m_fnew.at(i)(
                destination.index) =
                f_destination[i]  + feq[i];
    }


/*	fe_boundary_values.getData().m_fnew.at(destination.direction)(
					destination.index) =
            f_destination[destination.direction]; */

}

// Explicit instantiation
template class ThermalBounceBack<2> ;
template class ThermalBounceBack<3> ;

} /* namespace natrium */
