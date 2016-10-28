/*
 * KBCCentral.cpp
 *
 *  Created on: 29.03.2016
 *      Author: Dominik Wilde
 */

#include "KBCCentral.h"
#define EVALUATE_GAMMA // if defined, an evaluation over time of the stabilizer gamma will be carried out

namespace natrium {

KBCCentral::KBCCentral(double relaxationParameter, double dt,
		const boost::shared_ptr<Stencil> stencil) :
		MRT(relaxationParameter, dt, stencil), m_D(setMRTWeights()), m_S(
				setRelaxationRates()) {

}

KBCCentral::~KBCCentral() {
}

void KBCCentral::collideAll(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const {

	if (Stencil_D2Q9 == getStencil()->getStencilType()) {
		collideAllD2Q9(f, densities, velocities, locally_owned_dofs,
				inInitializationProcedure);
		/*else if (Stencil_D3Q19 == getStencil()->getStencilType()) {
		 collideAllD3Q19(f, densities, velocities, locally_owned_dofs,
		 inInitializationProcedure);*/

	} else {
		throw CollisionException("KBC_Central only implemented for D2Q9");
	}
}

void KBCCentral::collideAllD2Q9(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const {


		size_t Q = getQ();

		double scaling = getStencil()->getScaling();
		double cs2 = getStencil()->getSpeedOfSoundSquare()/scaling/scaling;
		vector<double> meq(Q); 	// moment equilibrium distribution functions
		vector<double> m(Q);   	// moment distribution

		//for all degrees of freedom on current processor
		dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
		dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());
		for (it = locally_owned_dofs.begin(); it != end; it++) {
			size_t i = *it;

		if (densities(i) < 1e-10) {
			throw CollisionException(
					"Densities too small (< 1e-10) for collisions. Decrease time step size.");
		}

			// transform the velocity space into moment space
			m.at(0) = f.at(0)(i) + f.at(1)(i) + f.at(2)(i) + f.at(3)(i)
					+ f.at(4)(i) + f.at(5)(i) + f.at(6)(i) + f.at(7)(i)
					+ f.at(8)(i);
			m.at(1) = -4 * f.at(0)(i) - f.at(1)(i) - f.at(2)(i) - f.at(3)(i)
					- f.at(4)(i)
					+ 2 * (f.at(5)(i) + f.at(6)(i) + f.at(7)(i) + f.at(8)(i));
//			m.at(2) = 4 * f.at(0)(i)
//					- 2 * (f.at(1)(i) + f.at(2)(i) + f.at(3)(i) + f.at(4)(i))
					+ f.at(5)(i) + f.at(6)(i) + f.at(7)(i) + f.at(8)(i);
			m.at(3) = f.at(1)(i) - f.at(3)(i) + f.at(5)(i) - f.at(6)(i)
					- f.at(7)(i) + f.at(8)(i);
//			m.at(4) = -2 * (f.at(1)(i) - f.at(3)(i)) + f.at(5)(i) - f.at(6)(i)
//					- f.at(7)(i) + f.at(8)(i);
			m.at(5) = f.at(2)(i) - f.at(4)(i) + f.at(5)(i) + f.at(6)(i)
					- f.at(7)(i) - f.at(8)(i);
//			m.at(6) = -2 * (f.at(2)(i) - f.at(4)(i)) + f.at(5)(i) + f.at(6)(i)
//					- f.at(7)(i) - f.at(8)(i);
			m.at(7) = f.at(1)(i) - f.at(2)(i) + f.at(3)(i) - f.at(4)(i);
			m.at(8) = f.at(5)(i) - f.at(6)(i) + f.at(7)(i) - f.at(8)(i);

			// calculate the moment equilibrium distribution function
			double rho = m.at(0);
			double jx = m.at(3);
			double jy = m.at(5);

			// calculate density
			densities(i) = rho;

			if (not inInitializationProcedure) {

				velocities.at(0)(i) = scaling / densities(i)
						* jx;
				velocities.at(1)(i) = scaling / densities(i)
						* jy;
			}




			meq.at(0) = rho;
			meq.at(1) = -2*rho+3/rho*(jx*jx+jy*jy);
			meq.at(2) = 0;
			meq.at(3) = jx;
			meq.at(4) = 0;
			meq.at(5) = jy;
			meq.at(6) = 0;
			meq.at(7) = jx*jx-jy*jy;
			meq.at(8) = jx*jy;


			//relax and rescale the moments

				m.at(1) = m.at(1) + -1./(getRelaxationParameter()+0.5) * (m.at(1) - meq.at(1));
				m.at(7) = m.at(7) + -1./(getRelaxationParameter()+0.5) * (m.at(7) - meq.at(7));
				m.at(8) = m.at(8) + -1./(getRelaxationParameter()+0.5) * (m.at(8) - meq.at(8));


//cout << getPrefactor() << " " << -1./(getRelaxationParameter()+0.5) << " " << getTime() << endl;

			m.at(2)=-rho-m.at(1);
			m.at(4)=-jx;
			m.at(6)=-jy;

			for (size_t j = 0; j < Q; j++) {
			m.at(j) /= m_D.at(j);
			}

			//transform the momentum space back into velocity space
			f.at(0)(i) = m.at(0) - 4 * (m.at(1) - m.at(2));
			f.at(1)(i) = m.at(0) - m.at(1) - 2 * (m.at(2) + m.at(4)) + m.at(3)
					+ m.at(7);
			f.at(2)(i) = m.at(0) - m.at(1) - 2 * (m.at(2) + m.at(6)) + m.at(5)
					- m.at(7);
			f.at(3)(i) = m.at(0) - m.at(1) - 2 * (m.at(2) - m.at(4)) - m.at(3)
					+ m.at(7);
			f.at(4)(i) = m.at(0) - m.at(1) - 2 * (m.at(2) - m.at(6)) - m.at(5)
					- m.at(7);
			f.at(5)(i) = m.at(0) + m.at(1) + m.at(1) + m.at(2) + m.at(3)
					+ m.at(4) + m.at(5) + m.at(6) + m.at(8);
			f.at(6)(i) = m.at(0) + m.at(1) + m.at(1) + m.at(2) - m.at(3)
					- m.at(4) + m.at(5) + m.at(6) - m.at(8);
			f.at(7)(i) = m.at(0) + m.at(1) + m.at(1) + m.at(2) - m.at(3)
					- m.at(4) - m.at(5) - m.at(6) + m.at(8);
			f.at(8)(i) = m.at(0) + m.at(1) + m.at(1) + m.at(2) + m.at(3)
					+ m.at(4) - m.at(5) - m.at(6) - m.at(8);

	}

//#ifdef EVALUATE_GAMMA
//	writeDeviation(gamma.getAverage(), gamma.getDeviation(),
//			entropy.getAverage(), entropy.getDeviation());
//#endif

}
}
