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
		MRT(relaxationParameter, dt, stencil), m_D(setMRTWeights()), m_S(
				setRelaxationRates()) {

	//cout << "test"<< endl;

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

		size_t n_dofs = f.at(0).size();
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

/*		if (densities(i) < 1e-10) {
			throw CollisionException(
					"Densities too small (< 1e-10) for collisions. Decrease time step size.");
		}*/

			if (not inInitializationProcedure) {

				velocities.at(0)(i) = scaling / densities(i)
						* (f.at(1)(i) + f.at(5)(i) + f.at(8)(i) - f.at(3)(i)
								- f.at(6)(i) - f.at(7)(i));
				velocities.at(1)(i) = scaling / densities(i)
						* (f.at(2)(i) + f.at(5)(i) + f.at(6)(i) - f.at(4)(i)
								- f.at(7)(i) - f.at(8)(i));
			}

			double ux = velocities.at(0)(i) / scaling;
			double uy = velocities.at(1)(i) / scaling;

			if (i == 1) {
			}

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



			for (int a = 0; a<9; a++)
			{
				//cout << "f" << a << " " << f.at(a)(i) << endl;
			}

			for (int a = 0; a<9; a++)
			{
				//cout << "m" << a << " " << m.at(a) << endl;
			}

			// calculate the moment equilibrium distribution function
			double rho = m.at(0);
			double e = m.at(1);
			double jx = m.at(3);
			double jy = m.at(5);
			//assert(jx-ux<1e-10);
			double eps = m.at(2);
			double qx = m.at(4);
			double qy = m.at(6);
			double pxx = m.at(7);
			double pxy = m.at(8);



			meq.at(0) = rho;
			meq.at(1) = jx*jy/rho;
			meq.at(2) = 0;
			meq.at(3) = jx;
			meq.at(4) = 0;
			meq.at(5) = jy;
			meq.at(6) = 0;
			meq.at(7) = 1./3.*rho+jx*jx/rho;
			meq.at(8) = 1./3.*rho+jy*jy/rho;

			for (int a = 0; a<9; a++)
			{
				//cout << "meq" << a << " " << meq.at(a) << endl;
			}

			//relax and rescale the moments
			for (size_t j = 0; j < Q; j++) {
				m.at(j) = m.at(j) + -0.0989903 * (m.at(j) - meq.at(j));
				//////cout << m.at(j);

			}
//cout << getPrefactor() << endl;

			m.at(2)=-rho-m.at(1);
			m.at(4)=-jx;
			m.at(6)=-jy;

			for (size_t j = 0; j < Q; j++) {
			m.at(j) /= m_D.at(j);
			}

			for (int a = 0; a<9; a++)
			{
				//cout << "mpc" << a << " " << m.at(a) << endl;
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

			for (int a = 0; a<9; a++)
			{
				//cout << "fpc" << a << " " << f.at(a)(i) << endl;
			}

		}
	}
} // namespace natrium

