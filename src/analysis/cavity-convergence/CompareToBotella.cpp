/*
 * CompareToBotella.cpp
 *
 *  Created on: 14.03.2016
 *      Author: akraem3m
 */

#include "CompareToBotella.h"

#include <fstream>

#include "deal.II/numerics/fe_field_function.h"
#include "deal.II/base/mpi.h"

namespace natrium {

CompareToBotella::CompareToBotella(const CFDSolver<2>& solver, size_t reynolds,
		string filename) :
		DataProcessor<2>(solver), m_Re(reynolds), m_y(17), m_referenceU(17), m_x(
				17), m_referenceV(17), m_fileName(filename) {
	m_x(0) = 1.0;
	m_x(1) = 0.9375;
	m_x(2) = 0.9297;
	m_x(3) = 0.9219;
	m_x(4) = 0.9062;
	m_x(5) = 0.8437;
	m_x(6) = 0.7734;
	m_x(7) = 0.7656;
	m_x(8) = 0.5000;
	m_x(9) = 0.1953;
	m_x(10) = 0.1406;
	m_x(11) = 0.0937;
	m_x(12) = 0.0547;
	m_x(13) = 0.0469;
	m_x(14) = 0.0391;
	m_x(15) = 0.0312;
	m_x(16) = 0.0;

	m_y(0) = 1.0000;
	m_y(1) = 0.9766;
	m_y(2) = 0.9688;
	m_y(3) = 0.9609;
	m_y(4) = 0.9531;
	m_y(5) = 0.8516;
	m_y(6) = 0.7344;
	m_y(7) = 0.6172;
	m_y(8) = 0.5000;
	m_y(9) = 0.4531;
	m_y(10) = 0.2813;
	m_y(11) = 0.1719;
	m_y(12) = 0.1016;
	m_y(13) = 0.0703;
	m_y(14) = 0.0625;
	m_y(15) = 0.0547;
	m_y(16) = 0.0000;
	makeReferenceU(m_referenceU, reynolds);
	makeReferenceV(m_referenceV, reynolds);

}

CompareToBotella::~CompareToBotella() {
	// TODO Auto-generated destructor stub
}

void CompareToBotella::apply() {
	if (m_solver.getIteration() % 1000 != 0)
		return;
	std::ofstream out(m_fileName);
	if (is_MPI_rank_0()) {
		out << "y   uy   uy_ref      x   vx  vx_ref" << endl;
	}

	m_uError = 0;
	m_vError = 0;

	// make fe field function to evaluate the solution
	dealii::Functions::FEFieldFunction<2, dealii::DoFHandler<2>,
			distributed_vector> fe_function_u(
			*m_solver.getAdvectionOperator()->getDoFHandler(),
			m_solver.getVelocity().at(0),
			m_solver.getAdvectionOperator()->getMapping());
	dealii::Functions::FEFieldFunction<2, dealii::DoFHandler<2>,
			distributed_vector> fe_function_v(
			*m_solver.getAdvectionOperator()->getDoFHandler(),
			m_solver.getVelocity().at(1),
			m_solver.getAdvectionOperator()->getMapping());

	// for all test points in Botella table
	for (size_t i = 0; i < 17; i++) {

		// generate test points
		dealii::Point<2> x_point(m_x(i), 0.5);
		dealii::Point<2> y_point(0.5, m_y(i));
		double ui;
		double vi;
		double ui_mpi;
		double vi_mpi;
		// evaluate solution and synchronize via MPI
		try {
			vi = fe_function_v.value(x_point);
		} catch (const dealii::VectorTools::ExcPointNotAvailableHere &) {
			vi = -1000000000;
		}
		vi_mpi = dealii::Utilities::MPI::max(vi, MPI_COMM_WORLD);

		try {
			ui = fe_function_u.value(y_point);
		} catch (const dealii::VectorTools::ExcPointNotAvailableHere &) {
			ui = -10000000000;
		}
		ui_mpi = dealii::Utilities::MPI::max(ui, MPI_COMM_WORLD);
		// calculate velocities at given points
		if (is_MPI_rank_0()) {
			out << m_y(i) << " " << ui_mpi << " " << m_referenceU(i) << " "
					<< m_x(i) << " " << vi_mpi << " " << m_referenceV(i)
					<< endl;
		}
		m_uError += (ui_mpi - m_referenceU(i)) * (ui_mpi - m_referenceU(i));
		m_vError += (vi_mpi - m_referenceV(i)) * (vi_mpi - m_referenceV(i));
	}
	out.close();

	m_uError = sqrt(m_uError / 17.);
	m_vError = sqrt(m_vError / 17.);
}

void CompareToBotella::printFinalVelocities() {
	std::stringstream final_name;
	final_name << m_fileName << "final";
	std::ofstream out(final_name.str());
	if (is_MPI_rank_0()) {
		out << "y   uy   x   vx " << endl;
	}

	// make fe field function to evaluate the solution
	dealii::Functions::FEFieldFunction<2, dealii::DoFHandler<2>,
			distributed_vector> fe_function_u(
			*m_solver.getAdvectionOperator()->getDoFHandler(),
			m_solver.getVelocity().at(0),
			m_solver.getAdvectionOperator()->getMapping());
	dealii::Functions::FEFieldFunction<2, dealii::DoFHandler<2>,
			distributed_vector> fe_function_v(
			*m_solver.getAdvectionOperator()->getDoFHandler(),
			m_solver.getVelocity().at(1),
			m_solver.getAdvectionOperator()->getMapping());

	// for all test points in Botella table
	for (size_t i = 0; i < 151; i++) {

		// generate test points
		dealii::Point<2> x_point(1.0 * i / 150., 0.5);
		dealii::Point<2> y_point(0.5, 1.0 * i / 150.);
		double ui;
		double vi;
		double ui_mpi;
		double vi_mpi;
		// evaluate solution and synchronize via MPI
		try {
			vi = fe_function_v.value(x_point);
		} catch (const dealii::VectorTools::ExcPointNotAvailableHere &) {
			vi = -1000000000;
		}
		vi_mpi = dealii::Utilities::MPI::max(vi, MPI_COMM_WORLD);

		try {
			ui = fe_function_u.value(y_point);
		} catch (const dealii::VectorTools::ExcPointNotAvailableHere &) {
			ui = -10000000000;
		}
		ui_mpi = dealii::Utilities::MPI::max(ui, MPI_COMM_WORLD);
		// calculate velocities at given points
		if (is_MPI_rank_0()) {
			out << y_point(1) << " " << ui_mpi << " " << x_point(0) << " "
					<< vi_mpi << endl;
		}
	}
	out.close();

}

void CompareToBotella::makeReferenceU(numeric_vector& u, size_t reynolds_number) {
	switch (reynolds_number) {
	// Hier waren die Vorzeichen umgedreht	
	case 1000: {
		u(0) = 1.0;
		u(1) = 0.6644227;
		u(2) = 0.5808359;
		u(3) = 0.5169277;
		u(4) = 0.4723329;
		u(5) = 0.3372212;
		u(6) = 0.1886747;
		u(7) = 0.0570178;
		u(8) = -0.0620561;
		u(9) = -0.1081999;
		u(10) = -0.2803696;
		u(11) = -0.3885691;
		u(12) = -0.3004561;
		u(13) = -0.2228955;
		u(14) = -0.2023300;
		u(15) = -0.1812881;
		u(16) = 0.0000000;
		return;
	}
	
	default: {
		return;
	}
	}
}

void CompareToBotella::makeReferenceV(numeric_vector& v, size_t reynolds_number) {
	switch (reynolds_number) {
	case 1000: {
		v(0) = 0.000;
		v(1) = -0.2279225;
		v(2) = -0.2936869;
		v(3) = -0.3553213;
		v(4) = -0.4103754;
		v(5) = -0.5264392;
		v(6) = -0.4264545;
		v(7) = -0.3202137;
		v(8) = 0.0257995;
		v(9) = 0.3253592;
		v(10) = 0.3339924;
		v(11) = 0.3769189;
		v(12) = 0.3330442;
		v(13) = 0.3099097;
		v(14) = 0.2962703;
		v(15) = 0.2807056;
		v(16) = 0.0000000;
		return;
	}
	default: {
		return;
	}
	}
}

} /* namespace natrium */
