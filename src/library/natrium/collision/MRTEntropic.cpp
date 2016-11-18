/*
 * MRTEntropic.cpp
 *
 *  Created on: 31.10.2016
 *      Author: Dominik Wilde
 */

#include "MRTEntropic.h"
#define EVALUATE_GAMMA // if defined, an evaluation over time of the stabilizer gamma will be carried out

namespace natrium {

MRTEntropic::MRTEntropic(double relaxationParameter, double dt,
		const boost::shared_ptr<Stencil> stencil) :
		MRT(relaxationParameter, dt, stencil) {

}

MRTEntropic::~MRTEntropic() {
}

void MRTEntropic::collideAll(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const {

	if (Stencil_D2Q9 == getStencil()->getStencilType()) {
		collideAllD2Q9(f, densities, velocities, locally_owned_dofs,
				inInitializationProcedure);
	} else if (Stencil_D3Q19 == getStencil()->getStencilType()) {
		collideAllD3Q19(f, densities, velocities, locally_owned_dofs,
				inInitializationProcedure);

	} else {
		throw CollisionException("KBC_Central only implemented for D2Q9/D3Q19");
	}
}

void MRTEntropic::collideAllD2Q9(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const {

	size_t Q = getQ();

	double scaling = getStencil()->getScaling();
	double cs2 = getStencil()->getSpeedOfSoundSquare() / scaling / scaling;
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
		m.at(0) = f.at(0)(i) + f.at(1)(i) + f.at(2)(i) + f.at(3)(i) + f.at(4)(i)
				+ f.at(5)(i) + f.at(6)(i) + f.at(7)(i) + f.at(8)(i);
		m.at(1) = -4 * f.at(0)(i) - f.at(1)(i) - f.at(2)(i) - f.at(3)(i)
				- f.at(4)(i)
				+ 2 * (f.at(5)(i) + f.at(6)(i) + f.at(7)(i) + f.at(8)(i));
//			m.at(2) = 4 * f.at(0)(i)
//					- 2 * (f.at(1)(i) + f.at(2)(i) + f.at(3)(i) + f.at(4)(i))
		+f.at(5)(i) + f.at(6)(i) + f.at(7)(i) + f.at(8)(i);
		m.at(3) = f.at(1)(i) - f.at(3)(i) + f.at(5)(i) - f.at(6)(i) - f.at(7)(i)
				+ f.at(8)(i);
//			m.at(4) = -2 * (f.at(1)(i) - f.at(3)(i)) + f.at(5)(i) - f.at(6)(i)
//					- f.at(7)(i) + f.at(8)(i);
		m.at(5) = f.at(2)(i) - f.at(4)(i) + f.at(5)(i) + f.at(6)(i) - f.at(7)(i)
				- f.at(8)(i);
//			m.at(6) = -2 * (f.at(2)(i) - f.at(4)(i)) + f.at(5)(i) + f.at(6)(i)
//					- f.at(7)(i) - f.at(8)(i);
		m.at(7) = f.at(1)(i) - f.at(2)(i) + f.at(3)(i) - f.at(4)(i);
		m.at(8) = f.at(5)(i) - f.at(6)(i) + f.at(7)(i) - f.at(8)(i);

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

		// calculate the moment equilibrium distribution function
		double rho = m.at(0);
		double jx = m.at(3);
		double jy = m.at(5);

		// calculate density
		densities(i) = rho;

		if (not inInitializationProcedure) {

			velocities.at(0)(i) = scaling / densities(i) * jx;
			velocities.at(1)(i) = scaling / densities(i) * jy;
		}

		meq.at(0) = rho;
		meq.at(1) = -2 * rho + 3 / rho * (jx * jx + jy * jy);
		meq.at(2) = 0;
		meq.at(3) = jx;
		meq.at(4) = 0;
		meq.at(5) = jy;
		meq.at(6) = 0;
		meq.at(7) = jx * jx - jy * jy;
		meq.at(8) = jx * jy;

		//relax and rescale the moments

		m.at(1) = m.at(1)
				+ -1. / (getRelaxationParameter() + 0.5)
						* (m.at(1) - meq.at(1));
		m.at(7) = m.at(7)
				+ -1. / (getRelaxationParameter() + 0.5)
						* (m.at(7) - meq.at(7));
		m.at(8) = m.at(8)
				+ -1. / (getRelaxationParameter() + 0.5)
						* (m.at(8) - meq.at(8));

//cout << getPrefactor() << " " << -1./(getRelaxationParameter()+0.5) << " " << getTime() << endl;

		m.at(2) = -rho - m.at(1);
		m.at(4) = -jx;
		m.at(6) = -jy;

		for (size_t j = 0; j < Q; j++) {
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
		f.at(5)(i) = m.at(0) + m.at(1) + m.at(1) + m.at(2) + m.at(3) + m.at(4)
				+ m.at(5) + m.at(6) + m.at(8);
		f.at(6)(i) = m.at(0) + m.at(1) + m.at(1) + m.at(2) - m.at(3) - m.at(4)
				+ m.at(5) + m.at(6) - m.at(8);
		f.at(7)(i) = m.at(0) + m.at(1) + m.at(1) + m.at(2) - m.at(3) - m.at(4)
				- m.at(5) - m.at(6) + m.at(8);
		f.at(8)(i) = m.at(0) + m.at(1) + m.at(1) + m.at(2) + m.at(3) + m.at(4)
				- m.at(5) - m.at(6) - m.at(8);

	}

//#ifdef EVALUATE_GAMMA
//	writeDeviation(gamma.getAverage(), gamma.getDeviation(),
//			entropy.getAverage(), entropy.getDeviation());
//#endif

}

void MRTEntropic::collideAllD3Q19(DistributionFunctions& f,
		distributed_vector& densities, vector<distributed_vector>& velocities,
		const dealii::IndexSet& locally_owned_dofs,
		bool inInitializationProcedure) const {

	double tm[19][19] = { { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
			1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 }, { -30.0, -11.0,
			-11.0, -11.0, -11.0, -11.0, -11.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0,
			8.0, 8.0, 8.0, 8.0, 8.0, 8.0 }, { 12.0, -4.0, -4.0, -4.0, -4.0,
			-4.0, -4.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
			1.0 }, { 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0,
			1.0, -1.0, 1.0, -1.0, 0, 0, 0, 0 }, { 0.0, -4.0, 4.0, 0.0, 0.0, 0.0,
			0.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 0, 0, 0, 0 }, {
			0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, 1.0, -1.0, -1.0, 0, 0, 0,
			0, 1.0, -1.0, 1.0, -1.0 }, { 0.0, 0.0, 0.0, -4.0, 4.0, 0.0, 0.0,
			1.0, 1.0, -1.0, -1.0, 0, 0, 0, 0, 1.0, -1.0, 1.0, -1.0 }, { 0.0,
			0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, -1.0,
			-1.0, 1.0, 1.0, -1.0, -1.0 }, { 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 4.0,
			0.0, 0.0, 0.0, 0.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0 }, {
			0.0, 2.0, 2.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
			1.0, 1.0, -2.0, -2.0, -2.0, -2.0 }, { 0.0, -4.0, -4.0, 2.0, 2.0,
			2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -2.0, -2.0, -2.0,
			-2.0 }, { 0.0, 0.0, 0.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0,
			-1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0, -2.0,
			-2.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 0.0,
			0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0,
			-1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			1.0, 0 - 1.0, -1.0, 1.0 }, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0 }, { 0.0,
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0,
			1.0, 0.0, 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			-1.0, -1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0 }, {
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
			-1.0, -1.0, -1.0, -1.0, 1.0, 1.0 } };

double invm[19][19]= {{1. / 19.,-5. / 399.,1. / 21.,0. / 1.,0. / 1.,0. / 1.,0. / 1.,0. / 1.,0. / 1.,0. / 1.,0. / 1.,0. / 1.,0. / 1.,0. / 1.,0. / 1.,0. / 1.,0. / 1.,0. / 1., 0. / 1.},
{1./ 19.,-11./ 2394.,-1./ 63.,1./ 10.,-1./ 10.,0./ 1.,0./ 1.,0./ 1.,0./ 1.,1./ 18.,-1./ 18.,0./ 1.,0./ 1.,0./ 1.,0./ 1.,0./ 1.,0./ 1.,0./ 1.,0./ 1.},
{1./ 19.,-11./ 2394.,-1./ 63.,-1./ 10.,1./ 10.,0./ 1.,0./ 1.,0./ 1.,0./ 1.,1./ 18.,-1./ 18.,0./ 1.,0./ 1.,0./ 1.,0./ 1.,0./ 1.,0./ 1.,0./ 1.,0./ 1.},
{1./ 19.,-11./ 2394.,-1./ 63.,0./ 1.,0./ 1.,1./ 10.,-1./ 10.,0./ 1.,0./ 1.,-1./ 36.,1./ 36.,1./ 12.,-1./ 12.,0./ 1.,0./ 1.,0./ 1.,0./ 1.,0./ 1.,0./ 1.},
{1./ 19.,-11./ 2394.,-1./ 63.,0./ 1.,0./ 1.,-1./ 10.,1./ 10.,0./ 1.,0./ 1.,-1./ 36.,1./ 36.,1./ 12.,-1./ 12.,0./ 1.,0./ 1.,0./ 1.,0./ 1.,0./ 1.,0./ 1.},
{1./ 19.,-11./ 2394.,-1./ 63.,0./ 1.,0./ 1.,0./ 1.,0./ 1.,1./ 10.,-1./ 10.,-1./ 36.,1./ 36.,-1./ 12.,1./ 12.,0./ 1.,0./ 1.,0./ 1.,0./ 1.,0./ 1.,0./ 1.},
{1./ 19.,-11./ 2394.,-1./ 63.,0./ 1.,0./ 1.,0./ 1.,0./ 1.,-1./ 10.,1./ 10.,-1./ 36.,1./ 36.,-1./ 12.,1./ 12.,0./ 1.,0./ 1.,0./ 1.,0./ 1.,0./ 1.,0./ 1.},
{1./ 19.,4./ 1197.,1./ 252.,1./ 10.,1./ 40.,1./ 10.,1./ 40.,0./ 1.,0./ 1.,1./ 36.,1./ 72.,1./ 12.,1./ 24.,1./ 4.,0./ 1.,0./ 1.,1./ 8.,-1./ 8.,0./ 1.},
{1./ 19.,4./ 1197.,1./ 252.,-1./ 10.,-1./ 40.,1./ 10.,1./ 40.,0./ 1.,0./ 1.,1./ 36.,1./ 72.,1./ 12.,1./ 24.,-1./ 4.,0./ 1.,0./ 1.,-1./ 8.,-1./ 8.,0./ 1.},
{1./ 19.,4./ 1197.,1./ 252.,1./ 10.,1./ 40.,-1./ 10.,-1./ 40.,0./ 1.,0./ 1.,1./ 36.,1./ 72.,1./ 12.,1./ 24.,-1./ 4.,0./ 1.,0./ 1.,1./ 8.,1./ 8.,0./ 1.},
{1./ 19.,4./ 1197.,1./ 252.,-1./ 10.,-1./ 40.,-1./ 10.,-1./ 40.,0./ 1.,0./ 1.,1./ 36.,1./ 72.,1./ 12.,1./ 24.,1./ 4.,0./ 1.,0./ 1.,-1./ 8.,1./ 8.,0./ 1.},
{1./ 19.,4./ 1197.,1./ 252.,1./ 10.,1./ 40.,0./ 1.,0./ 1.,1./ 10.,1./ 40.,1./ 36.,1./ 72.,-1./ 12.,-1./ 24.,0./ 1.,0./ 1.,1./ 4.,-1./ 8.,0./ 1.,1./ 8.},
{1./ 19.,4./ 1197.,1./ 252.,-1./ 10.,-1./ 40.,0./ 1.,0./ 1.,1./ 10.,1./ 40.,1./ 36.,1./ 72.,-1./ 12.,-1./ 24.,0./ 1.,0./ 1.,-1./ 4.,1./ 8.,0./ 1.,1./ 8.},
{1./ 19.,4./ 1197.,1./ 252.,1./ 10.,1./ 40.,0./ 1.,0./ 1.,-1./ 10.,-1./ 40.,1./ 36.,1./ 72.,-1./ 12.,-1./ 24.,0./ 1.,0./ 1.,-1./ 4.,-1./ 8.,0./ 1.,-1./ 8.},
{1./ 19.,4./ 1197.,1./ 252.,-1./ 10.,-1./ 40.,0./ 1.,0./ 1.,-1./ 10.,-1./ 40.,1./ 36.,1./ 72.,-1./ 12.,-1./ 24.,0./ 1.,0./ 1.,1./ 4.,1./ 8.,0./ 1.,-1./ 8.},
{1./ 19.,4./ 1197.,1./ 252.,0./ 1.,0./ 1.,1./ 10.,1./ 40.,1./ 10.,1./ 40.,-1./ 18.,-1./ 36.,0./ 1.,0./ 1.,0./ 1.,1./ 4.,0./ 1.,0./ 1.,1./ 8.,-1./ 8.},
{1./ 19.,4./ 1197.,1./ 252.,0./ 1.,0./ 1.,-1./ 10.,-1./ 40.,1./ 10.,1./ 40.,-1./ 18.,-1./ 36.,0./ 1.,0./ 1.,0./ 1.,-1./ 4.,0./ 1.,0./ 1.,-1./ 8.,-1./ 8.},
{1./ 19.,4./ 1197.,1./ 252.,0./ 1.,0./ 1.,1./ 10.,1./ 40.,-1./ 10.,-1./ 40.,-1./ 18.,-1./ 36.,0./ 1.,0./ 1.,0./ 1.,-1./ 4.,0./ 1.,0./ 1.,1./ 8.,1./ 8.},
{1./ 19.,4./ 1197.,1./ 252.,0./ 1.,0./ 1.,-1./ 10.,-1./ 40.,-1./ 10.,-1./ 40.,-1./ 18.,-1./ 36.,0./ 1.,0./ 1.,0./ 1.,1./ 4.,0./ 1.,0./ 1.,-1./ 8.,1./ 8.}};

size_t Q = getQ();

double scaling = getStencil()->getScaling();
double cs2 = getStencil()->getSpeedOfSoundSquare() / scaling / scaling;
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

	for (int p = 0; p < 19; p++) {
		m[p]=0;
		for (int q = 0; q < 19; q++) {
			m[p] += tm[p][q] * f.at(q)(i);
		}
	}

	double rho = m.at(0);
	double jx = m.at(3);
	double jy = m.at(5);
	double jz = m.at(7);

	densities(i) = rho;

	if (not inInitializationProcedure) {

		velocities.at(0)(i) = scaling / densities(i) * jx;
		velocities.at(1)(i) = scaling / densities(i) * jy;
		velocities.at(2)(i) = scaling / densities(i) * jz;
	}

	meq.at(1) = -11 * rho + 19. / rho * (jx * jx + jy * jy + jz * jz);
	meq.at(9) = 1. / rho * (2 * jx * jx - (jy * jy + jz * jz));
	meq.at(11) = 1. / rho * (jy * jy - jz * jz);
	meq.at(13) = 1. / rho * jx * jy;
	meq.at(14) = 1. / rho * jy * jz;
	meq.at(15) = 1. / rho * jx * jz;

	m.at(1) = m.at(1)
	+ -1. / (getRelaxationParameter() + 0.5)
	* (m.at(1) - meq.at(1));
	m.at(9) = m.at(9)
	+ -1. / (getRelaxationParameter() + 0.5)
	* (m.at(9) - meq.at(9));
	m.at(11) = m.at(11)
	+ -1. / (getRelaxationParameter() + 0.5)
	* (m.at(11) - meq.at(11));
	m.at(13) = m.at(13)
	+ -1. / (getRelaxationParameter() + 0.5)
	* (m.at(13) - meq.at(13));
	m.at(14) = m.at(14)
	+ -1. / (getRelaxationParameter() + 0.5)
	* (m.at(14) - meq.at(14));
	m.at(15) = m.at(15)
	+ -1. / (getRelaxationParameter() + 0.5)
	* (m.at(15) - meq.at(15));

	m.at(2) = -7. / 38 * rho - 11. / 38 * m.at(1);
	m.at(4) = -2. / 3. * jx;
	m.at(6) = -2. / 3. * jy;
	m.at(8) = -2. / 3. * jz;
	m.at(10) = -1. / 2. * m.at(9);
	m.at(12) = -1. / 2. * m.at(11);
	m.at(16) = 0;
	m.at(17) = 0;
	m.at(18) = 0;

	for (int p = 0; p < 19; p++) {
		f.at(p)(i)=0;
		for (int q = 0; q < 19; q++) {
			f.at(p)(i) += invm[p][q] * m.at(q);
		}
	}

}
}
}
