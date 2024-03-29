/*
 * MRTStandard.cpp
 *
 *  Created on: 25.07.2015
 *      Author: dominik
 */

#include "MRTStandard.h"

namespace natrium {

MRTStandard::MRTStandard(double relaxationParameter, double dt,
		const boost::shared_ptr<Stencil> stencil) :
		MRT(relaxationParameter, dt, stencil) {

}

MRTStandard::~MRTStandard() {
}

void MRTStandard::collideAll(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const {

	if (Stencil_D2Q9 == getStencil()->getStencilType()) {
		collideAllD2Q9(f, densities, velocities, locally_owned_dofs,
				inInitializationProcedure);
	} /*else if (Stencil_D3Q19 == getStencil()->getStencilType()) {
		collideAllD3Q19(f, densities, velocities, locally_owned_dofs,
				inInitializationProcedure);
	} else if (Stencil_D3Q15 == getStencil()->getStencilType()) {
		collideAllD3Q15(f, densities, velocities, locally_owned_dofs,
				inInitializationProcedure);
	} */else {
		throw CollisionException("MRT only implemented for D2Q9");
		// Inefficient collision
		//BGK::collideAll(f, densities, velocities, locally_owned_dofs,
		//		inInitializationProcedure);
	}
}

void MRTStandard::collideAllD2Q9(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const {

		//size_t n_dofs = f.at(0).size();
		size_t Q = getQ();
		double scaling = getStencil()->getScaling();
		vector<double> meq(Q); 	// moment equilibrium distribution functions
		vector<double> m(Q);   	// moment distribution

		//for all degrees of freedom on current processor
		dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
		dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());
		for (it = locally_owned_dofs.begin(); it != end; it++) {
			size_t i = *it;

			// calculate density
			densities(i) = f.at(0)(i) + f.at(1)(i) + f.at(2)(i) + f.at(3)(i)
					+ f.at(4)(i) + f.at(5)(i) + f.at(6)(i) + f.at(7)(i)
					+ f.at(8)(i);

			if (not inInitializationProcedure) {

				velocities.at(0)(i) = scaling / densities(i)
						* (f.at(1)(i) + f.at(5)(i) + f.at(8)(i) - f.at(3)(i)
								- f.at(6)(i) - f.at(7)(i));
				velocities.at(1)(i) = scaling / densities(i)
						* (f.at(2)(i) + f.at(5)(i) + f.at(6)(i) - f.at(4)(i)
								- f.at(7)(i) - f.at(8)(i));
			}


//cout<< "Test"<< endl;
			vector<double> MRTWeights(9);
			MRTWeights.at(0) = 9.0;
			MRTWeights.at(1) = 36.0;
			MRTWeights.at(2) = 36.0;
			MRTWeights.at(3) = 6.0;
			MRTWeights.at(4) = 12.0;
			MRTWeights.at(5) = 6.0;
			MRTWeights.at(6) = 12.0;
			MRTWeights.at(7) = 4.0;
			MRTWeights.at(8) = 4.0;


				vector<double> s(9);
				s.at(0) = s.at(3) = s.at(5) = 0.0;
				s.at(7) = s.at(8) =  1. / (getRelaxationParameter() + 0.5);
				s.at(4) = s.at(6) = 8.0 * (2.0 - s.at(7)) / (8.0 - s.at(7));

				s.at(1) = 1.6;
				s.at(2) = 1.8;


			// transform the velocity space into moment space
			m.at(0) = f.at(0)(i) + f.at(1)(i) + f.at(2)(i) + f.at(3)(i)
					+ f.at(4)(i) + f.at(5)(i) + f.at(6)(i) + f.at(7)(i)
					+ f.at(8)(i);
			m.at(1) = -4 * f.at(0)(i) - f.at(1)(i) - f.at(2)(i) - f.at(3)(i)
					- f.at(4)(i)
					+ 2 * (f.at(5)(i) + f.at(6)(i) + f.at(7)(i) + f.at(8)(i));
			m.at(2) = 4 * f.at(0)(i)
					- 2 * (f.at(1)(i) + f.at(2)(i) + f.at(3)(i) + f.at(4)(i))
					+ f.at(5)(i) + f.at(6)(i) + f.at(7)(i) + f.at(8)(i);
			m.at(3) = f.at(1)(i) - f.at(3)(i) + f.at(5)(i) - f.at(6)(i)
					- f.at(7)(i) + f.at(8)(i);
			m.at(4) = -2 * (f.at(1)(i) - f.at(3)(i)) + f.at(5)(i) - f.at(6)(i)
					- f.at(7)(i) + f.at(8)(i);
			m.at(5) = f.at(2)(i) - f.at(4)(i) + f.at(5)(i) + f.at(6)(i)
					- f.at(7)(i) - f.at(8)(i);
			m.at(6) = -2 * (f.at(2)(i) - f.at(4)(i)) + f.at(5)(i) + f.at(6)(i)
					- f.at(7)(i) - f.at(8)(i);
			m.at(7) = f.at(1)(i) - f.at(2)(i) + f.at(3)(i) - f.at(4)(i);
			m.at(8) = f.at(5)(i) - f.at(6)(i) + f.at(7)(i) - f.at(8)(i);

			double ux = m.at(3);
			double uy = m.at(5);
			//cout<< "Test2"<< endl;
			// calculate the moment equilibrium distribution function

			meq.at(0) = densities(i);
			meq.at(1) = densities(i) * (-2 + 3 * (ux * ux + uy * uy));
			meq.at(2) = densities(i) * (1 - 3 * (ux * ux + uy * uy));
			meq.at(3) = densities(i) * ux;
			meq.at(4) = -densities(i) * ux;
			meq.at(5) = densities(i) * uy;
			meq.at(6) = -densities(i) * uy;
			meq.at(7) = densities(i) * (ux * ux - uy * uy);
			meq.at(8) = densities(i) * ux * uy;

			//relax and rescale the moments
			for (size_t j = 0; j < 9; j++) {
				m.at(j) = m.at(j) - s.at(j) * (m.at(j) - meq.at(j));
				m.at(j) /= MRTWeights.at(j);
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

			//cout<< "Test3"<< endl;

		}
	}
} // namespace natrium

