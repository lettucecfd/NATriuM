/*
 * BGKMultistep.cpp
 *
 *  Created on: 20.06.2016
 *      Author: dominik
 */

#include "BGKMultistep.h"
#include "../solver/DistributionFunctions.h"

namespace natrium {

BGKMultistep::BGKMultistep(double relaxationParameter, double dt,
		const boost::shared_ptr<Stencil> stencil, int model) :
		BGK(relaxationParameter, dt, stencil), MultistepCollisionData() {
	// TODO Auto-generated constructor stub
	setTimeStep(dt);
	if (model==0)
	{m_model = ADAMSMOULTON4;
	cout << "AM4 selected";}
	if (model==1)
	{m_model = BDF2;
	cout << "BDF2 selected";}

}

BGKMultistep::~BGKMultistep() {
	// TODO Auto-generated destructor stub
}

double BGKMultistep::getEquilibriumDistribution(size_t i,
		const numeric_vector& u, const double rho) const {

	assert(i < getStencil()->getQ());
	assert(rho > 0);
	assert(u.size() == getStencil()->getD());
	assert(u(0) < 1000000000000000.);
	assert(u(1) < 1000000000000000.);

	double prefactor = getStencil()->getWeight(i) * rho;
	double uSquareTerm = -(u * u) / (2 * getStencil()->getSpeedOfSoundSquare());
	if (0 == i) {
		return prefactor * (1 + uSquareTerm);
	}
	double mixedTerm = (u * getStencil()->getDirection(i))
			/ getStencil()->getSpeedOfSoundSquare();
	return prefactor * (1 + mixedTerm * (1 + 0.5 * (mixedTerm)) + uSquareTerm);

} /// getEquilibriumDistribution

void BGKMultistep::collideAll(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const {

	if (m_firstCollision) {
		//cout << " First collision check" << endl ;
		m_formerF.reinit(getStencil()->getQ(), locally_owned_dofs,
		MPI_COMM_WORLD);
		m_formerFEq.reinit(getStencil()->getQ(), locally_owned_dofs,
		MPI_COMM_WORLD);
	}

	/*if (m_formerF.size() == 0){
	 m_formerF.reinit(m_stencil->getQ(), locally_owned_dofs , MPI_COMM_WORLD);
	 }

	 if (m_formerFEq.size() == 0){
	 m_formerFEq.reinit(m_stencil->getQ(), locally_owned_dofs , MPI_COMM_WORLD);
	 } */

	if (Stencil_D2Q9 == getStencil()->getStencilType()) {
		collideAllD2Q9(f, densities, velocities, locally_owned_dofs,
				inInitializationProcedure);
	}

	m_firstCollision = 0;
}

void BGKMultistep::collideAllD2Q9(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const {

	double cs2 = getStencil()->getSpeedOfSoundSquare(); ///(scaling * scaling);
	double lambda = getViscosity() / cs2;
	double multistep_factor0 =0;
	double multistep_factor1 =0;

	if (m_model == ADAMSMOULTON4) {
		multistep_factor0 = -1.
				/ (12. / 13. * lambda / getTimeStep() + 5. / 13.);
		multistep_factor1 = +1. / (12 * lambda / getTimeStep() + 5);
	}

	if (m_model == BDF2) {
		multistep_factor0 = -1. / (9. / 8. * lambda / getTimeStep() + 3. / 4.);
		multistep_factor1 = +1. / (9. / 2. * lambda / getTimeStep() + 3);
	}

	// Efficient collision for D2Q9
	size_t Q = 9;
	size_t D = 2;
	double scaling = getStencil()->getScaling();
	double prefactor = scaling / cs2;
	double relax_factor = -1 / (lambda / getTimeStep() + 0.5);

	// For the first collision, a BGK single step must be executed
	if (m_firstCollision) {
		multistep_factor0 = relax_factor;
		multistep_factor1 = 0;
	}

	//cout << "multi relax_factor: " << relax_factor << " prefactor: " << prefactor << endl;

	double dt = getTimeStep();

	// External force information
	ForceType force_type = getForceType();
	double force_x = getForceX();
	double force_y = getForceY();

	assert(f.size() == Q);
	assert(velocities.size() == D);

#ifdef DEBUG
	size_t n_dofs = f.at(0).size();
	for (size_t i = 0; i < Q; i++) {
		assert (f.at(i).size() == n_dofs);
	}
	assert (densities.size() == n_dofs);
	for (size_t i = 0; i < D; i++) {
		assert (velocities.at(i).size() == n_dofs);
	}
#endif

	// allocation
	double feq[9] = { };
	double scalar_product;
	double uSquareTerm;
	double mixedTerm;
	double weighting;
	double rho_i;
	double u_0_i;
	double u_1_i;
	double f_i[9] = { };
	double formerF_i[9] = { };
	double formerFEq_i[9] = { };

	distributed_vector& f0 = f.at(0);
	distributed_vector& f1 = f.at(1);
	distributed_vector& f2 = f.at(2);
	distributed_vector& f3 = f.at(3);
	distributed_vector& f4 = f.at(4);
	distributed_vector& f5 = f.at(5);
	distributed_vector& f6 = f.at(6);
	distributed_vector& f7 = f.at(7);
	distributed_vector& f8 = f.at(8);

	distributed_vector& former0 = m_formerF.at(0);
	distributed_vector& former1 = m_formerF.at(1);
	distributed_vector& former2 = m_formerF.at(2);
	distributed_vector& former3 = m_formerF.at(3);
	distributed_vector& former4 = m_formerF.at(4);
	distributed_vector& former5 = m_formerF.at(5);
	distributed_vector& former6 = m_formerF.at(6);
	distributed_vector& former7 = m_formerF.at(7);
	distributed_vector& former8 = m_formerF.at(8);

	distributed_vector& formerFEq0 = m_formerFEq.at(0);
	distributed_vector& formerFEq1 = m_formerFEq.at(1);
	distributed_vector& formerFEq2 = m_formerFEq.at(2);
	distributed_vector& formerFEq3 = m_formerFEq.at(3);
	distributed_vector& formerFEq4 = m_formerFEq.at(4);
	distributed_vector& formerFEq5 = m_formerFEq.at(5);
	distributed_vector& formerFEq6 = m_formerFEq.at(6);
	distributed_vector& formerFEq7 = m_formerFEq.at(7);
	distributed_vector& formerFEq8 = m_formerFEq.at(8);

	//for all degrees of freedom on current processor
	dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
	dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());
	for (it = locally_owned_dofs.begin(); it != end; it++) {
		size_t i = *it;

		// copy fi
		f_i[0] = f0(i);
		f_i[1] = f1(i);
		f_i[2] = f2(i);
		f_i[3] = f3(i);
		f_i[4] = f4(i);
		f_i[5] = f5(i);
		f_i[6] = f6(i);
		f_i[7] = f7(i);
		f_i[8] = f8(i);

		formerF_i[0] = former0(i);
		formerF_i[1] = former1(i);
		formerF_i[2] = former2(i);
		formerF_i[3] = former3(i);
		formerF_i[4] = former4(i);
		formerF_i[5] = former5(i);
		formerF_i[6] = former6(i);
		formerF_i[7] = former7(i);
		formerF_i[8] = former8(i);

		formerFEq_i[0] = formerFEq0(i);
		formerFEq_i[1] = formerFEq1(i);
		formerFEq_i[2] = formerFEq2(i);
		formerFEq_i[3] = formerFEq3(i);
		formerFEq_i[4] = formerFEq4(i);
		formerFEq_i[5] = formerFEq5(i);
		formerFEq_i[6] = formerFEq6(i);
		formerFEq_i[7] = formerFEq7(i);
		formerFEq_i[8] = formerFEq8(i);

		// calculate density
		rho_i = f_i[0] + f_i[1] + f_i[2] + f_i[3] + f_i[4] + f_i[5] + f_i[6]
				+ f_i[7] + f_i[8];
		densities(i) = rho_i;
		if (rho_i < 1e-10) {
			throw CollisionException(
					"Densities too small (< 1e-10) for collisions. Decrease time step size.");
		}

		if (not inInitializationProcedure) {
			// calculate macroscopic velocity (velocities.at()(i)) and equilibrium velocity (u_0_i, u_1_i)
			// for all velocity components
			u_0_i = scaling / rho_i
					* (f_i[1] + f_i[5] + f_i[8] - f_i[3] - f_i[6] - f_i[7]);
			u_1_i = scaling / rho_i
					* (f_i[2] + f_i[5] + f_i[6] - f_i[4] - f_i[7] - f_i[8]);

			if (force_type == NO_FORCING) {
				velocities.at(0)(i) = u_0_i;
				velocities.at(1)(i) = u_1_i;
			} else if (force_type == SHIFTING_VELOCITY) {
				velocities.at(0)(i) = u_0_i + 0.5 * dt * force_x / rho_i;
				velocities.at(1)(i) = u_1_i + 0.5 * dt * force_y / rho_i;
				u_0_i -= (double) 1 / relax_factor * dt * force_x / rho_i;
				u_1_i -= (double) 1 / relax_factor * dt * force_y / rho_i;
			} else if (force_type == EXACT_DIFFERENCE) {
				velocities.at(0)(i) = u_0_i + 0.5 * dt * force_x / rho_i;
				velocities.at(1)(i) = u_1_i + 0.5 * dt * force_y / rho_i;
			} else { // GUO
				u_0_i += 0.5 * dt * force_x / rho_i;
				u_1_i += 0.5 * dt * force_y / rho_i;
				velocities.at(0)(i) = u_0_i;
				velocities.at(1)(i) = u_1_i;
			}
		}

		// calculate equilibrium distribution
		scalar_product = u_0_i * u_0_i + u_1_i * u_1_i;
		uSquareTerm = -scalar_product / (2 * cs2);
		// direction 0
		weighting = 4. / 9. * rho_i;
		feq[0] = weighting * (1 + uSquareTerm);
		// directions 1-4
		weighting = 1. / 9. * rho_i;
		mixedTerm = prefactor * (u_0_i);
		feq[1] = weighting
				* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
		feq[3] = weighting
				* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
		mixedTerm = prefactor * (u_1_i);
		feq[2] = weighting
				* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
		feq[4] = weighting
				* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
		// directions 5-8
		weighting = 1. / 36. * rho_i;
		mixedTerm = prefactor * (u_0_i + u_1_i);
		feq[5] = weighting
				* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
		feq[7] = weighting
				* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
		mixedTerm = prefactor * (-u_0_i + u_1_i);
		feq[6] = weighting
				* (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
		feq[8] = weighting
				* (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);

		// BGK collision
		if (m_model == ADAMSMOULTON4) {
			f_i[0] += multistep_factor0 * (f_i[0] - feq[0])
					+ multistep_factor1 * (formerF_i[0] - formerFEq_i[0]);
			f_i[1] += multistep_factor0 * (f_i[1] - feq[1])
					+ multistep_factor1 * (formerF_i[1] - formerFEq_i[1]);
			f_i[2] += multistep_factor0 * (f_i[2] - feq[2])
					+ multistep_factor1 * (formerF_i[2] - formerFEq_i[2]);
			f_i[3] += multistep_factor0 * (f_i[3] - feq[3])
					+ multistep_factor1 * (formerF_i[3] - formerFEq_i[3]);
			f_i[4] += multistep_factor0 * (f_i[4] - feq[4])
					+ multistep_factor1 * (formerF_i[4] - formerFEq_i[4]);
			f_i[5] += multistep_factor0 * (f_i[5] - feq[5])
					+ multistep_factor1 * (formerF_i[5] - formerFEq_i[5]);
			f_i[6] += multistep_factor0 * (f_i[6] - feq[6])
					+ multistep_factor1 * (formerF_i[6] - formerFEq_i[6]);
			f_i[7] += multistep_factor0 * (f_i[7] - feq[7])
					+ multistep_factor1 * (formerF_i[7] - formerFEq_i[7]);
			f_i[8] += multistep_factor0 * (f_i[8] - feq[8])
					+ multistep_factor1 * (formerF_i[8] - formerFEq_i[8]);
		}

		if (m_model == BDF2 && !m_firstCollision) {
			f_i[0] =  4./3.*f_i[0] - 1./3. * formerF_i[0]+ multistep_factor0 * (f_i[0] - feq[0]) +  multistep_factor1 * (formerF_i[0] - formerFEq_i[0]);
			f_i[1] =  4./3.*f_i[1] - 1./3. * formerF_i[1]+ multistep_factor0 * (f_i[1] - feq[1]) +  multistep_factor1 * (formerF_i[1] - formerFEq_i[1]);
			f_i[2] =  4./3.*f_i[2] - 1./3. * formerF_i[2]+ multistep_factor0 * (f_i[2] - feq[2]) +  multistep_factor1 * (formerF_i[2] - formerFEq_i[2]);
			f_i[3] =  4./3.*f_i[3] - 1./3. * formerF_i[3]+ multistep_factor0 * (f_i[3] - feq[3]) +  multistep_factor1 * (formerF_i[3] - formerFEq_i[3]);
			f_i[4] =  4./3.*f_i[4] - 1./3. * formerF_i[4]+ multistep_factor0 * (f_i[4] - feq[4]) +  multistep_factor1 * (formerF_i[4] - formerFEq_i[4]);
			f_i[5] =  4./3.*f_i[5] - 1./3. * formerF_i[5]+ multistep_factor0 * (f_i[5] - feq[5]) +  multistep_factor1 * (formerF_i[5] - formerFEq_i[5]);
			f_i[6] =  4./3.*f_i[6] - 1./3. * formerF_i[6]+ multistep_factor0 * (f_i[6] - feq[6]) +  multistep_factor1 * (formerF_i[6] - formerFEq_i[6]);
			f_i[7] =  4./3.*f_i[7] - 1./3. * formerF_i[7]+ multistep_factor0 * (f_i[7] - feq[7]) +  multistep_factor1 * (formerF_i[7] - formerFEq_i[7]);
			f_i[8] =  4./3.*f_i[8] - 1./3. * formerF_i[8]+ multistep_factor0 * (f_i[8] - feq[8]) +  multistep_factor1 * (formerF_i[8] - formerFEq_i[8]);
		}

		if (m_model == BDF2 && m_firstCollision)
		{
			f_i[0] += relax_factor * (f_i[0] - feq[0]);
			f_i[1] += relax_factor * (f_i[1] - feq[1]);
			f_i[2] += relax_factor * (f_i[2] - feq[2]);
			f_i[3] += relax_factor * (f_i[3] - feq[3]);
			f_i[4] += relax_factor * (f_i[4] - feq[4]);
			f_i[5] += relax_factor * (f_i[5] - feq[5]);
			f_i[6] += relax_factor * (f_i[6] - feq[6]);
			f_i[7] += relax_factor * (f_i[7] - feq[7]);
			f_i[8] += relax_factor * (f_i[8] - feq[8]);
		}

		//	cout << "multistep_factor0" <<  multistep_factor0 << " ;  multistep_factor1: " << multistep_factor1 << endl;

		// Add Source term
		// Exact difference method (Kupershtokh)
		if (force_type == EXACT_DIFFERENCE) {
			ExternalForceFunctions::applyExactDifferenceForcingD2Q9(f_i,
					force_x, force_y, u_0_i, u_1_i, rho_i, getDt(), prefactor);
		} else if (force_type == GUO) {
			ExternalForceFunctions::applyGuoForcingD2Q9(f_i, force_x, force_y,
					u_0_i, u_1_i, -relax_factor, prefactor, dt);
		}

		//copy the former PDF
		former0(i) = f0(i);
		former1(i) = f1(i);
		former2(i) = f2(i);
		former3(i) = f3(i);
		former4(i) = f4(i);
		former5(i) = f5(i);
		former6(i) = f6(i);
		former7(i) = f7(i);
		former8(i) = f8(i);

		// copy the Equilibrium PDF
		formerFEq0(i) = feq[0];
		formerFEq1(i) = feq[1];
		formerFEq2(i) = feq[2];
		formerFEq3(i) = feq[3];
		formerFEq4(i) = feq[4];
		formerFEq5(i) = feq[5];
		formerFEq6(i) = feq[6];
		formerFEq7(i) = feq[7];
		formerFEq8(i) = feq[8];

		// copy to global variables
		f0(i) = f_i[0];
		f1(i) = f_i[1];
		f2(i) = f_i[2];
		f3(i) = f_i[3];
		f4(i) = f_i[4];
		f5(i) = f_i[5];
		f6(i) = f_i[6];
		f7(i) = f_i[7];
		f8(i) = f_i[8];

	} /* for all dofs */

} /* collideAllD2Q9 */

}/* namespace natrium */
