/*
 * CompareToGhia.cpp
 *
 *  Created on: 14.03.2016
 *      Author: akraem3m
 */

#include "CompareToGhia.h"

#include <fstream>

#include "deal.II/numerics/fe_field_function.h"
#include "deal.II/base/mpi.h"

namespace natrium {

CompareToGhia::CompareToGhia(const CFDSolver<2>& solver, size_t reynolds,
		string filename) :
		DataProcessor<2>(solver), m_Re(reynolds), m_y(17), m_referenceU(17), m_x(
				17), m_referenceV(17), m_fileName(filename) {
	m_x(0) = 1.0;
	m_x(1) = 0.9688;
	m_x(2) = 0.9609;
	m_x(3) = 0.9531;
	m_x(4) = 0.9453;
	m_x(5) = 0.9063;
	m_x(6) = 0.8594;
	m_x(7) = 0.8047;
	m_x(8) = 0.5000;
	m_x(9) = 0.2344;
	m_x(10) = 0.2266;
	m_x(11) = 0.1563;
	m_x(12) = 0.0938;
	m_x(13) = 0.0781;
	m_x(14) = 0.0703;
	m_x(15) = 0.0625;
	m_x(16) = 0.0000;
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

CompareToGhia::~CompareToGhia() {
	// TODO Auto-generated destructor stub
}

void CompareToGhia::apply() {
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

	// for all test points in ghia table
	for (size_t i = 0; i < 17; i++) {

		// generate test points
		dealii::Point<2> x_point(m_x(i), 0.5);
		dealii::Point<2> y_point(0.5, m_y(i));
		double ui = -1000000000;
		double vi = -1000000000;
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

void CompareToGhia::printFinalVelocities() {
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

	// for all test points in ghia table
	for (size_t i = 0; i < 151; i++) {

		// generate test points
		dealii::Point<2> x_point(1.0 * i / 150., 0.5);
		dealii::Point<2> y_point(0.5, 1.0 * i / 150.);
		double ui = -1000000000;
		double vi = -1000000000;
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

void CompareToGhia::makeReferenceU(numeric_vector& u, size_t reynolds_number) {
	switch (reynolds_number) {
	case 100: {
		u(0) = 1.00000;
		u(1) = 0.84123;
		u(2) = 0.78871;
		u(3) = 0.73722;
		u(4) = 0.68717;
		u(5) = 0.23151;
		u(6) = 0.00332;
		u(7) = -0.13641;
		u(8) = -0.20581;
		u(9) = -0.21090;
		u(10) = -0.15662;
		u(11) = -0.10150;
		u(12) = -0.06434;
		u(13) = -0.04775;
		u(14) = -0.04192;
		u(15) = -0.03717;
		u(16) = 0.00000;
		return;
	}
	case 400: {
		u(0) = 1.00000;
		u(1) = 0.75837;
		u(2) = 0.68439;
		u(3) = 0.61756;
		u(4) = 0.55892;
		u(5) = 0.29093;
		u(6) = 0.16256;
		u(7) = 0.02135;
		u(8) = -0.11477;
		u(9) = -0.17119;
		u(10) = -0.32726;
		u(11) = -0.24299;
		u(12) = -0.14612;
		u(13) = -0.10338;
		u(14) = -0.09266;
		u(15) = -0.08186;
		u(16) = 0.;
		return;
	}
	case 1000: {
		u(0) = 1.0;
		u(1) = 0.65928;
		u(2) = 0.57492;
		u(3) = 0.51117;
		u(4) = 0.46604;
		u(5) = 0.33304;
		u(6) = 0.18719;
		u(7) = 0.05702;
		u(8) = -0.06080;
		u(9) = -0.10648;
		u(10) = -0.27805;
		u(11) = -0.38289;
		u(12) = -0.29730;
		u(13) = -0.22220;
		u(14) = -0.20196;
		u(15) = -0.18109;
		u(16) = 0.00000;
		return;
	}
	case 3200: {
		u(0) = 1.0;
		u(1) = 0.53236;
		u(2) = 0.48296;
		u(3) = 0.46547;
		u(4) = 0.46101;
		u(5) = 0.34682;
		u(6) = 0.19791;
		u(7) = 0.07156;
		u(8) = -0.04272;
		u(9) = -0.86636;
		u(10) = -0.24427;
		u(11) = -0.34323;
		u(12) = -0.41933;
		u(13) = -0.37827;
		u(14) = -0.35344;
		u(15) = -0.32407;
		u(16) = 0.0000;
		return;
	}
	case 5000: {
		u(0) = 1.0;
		u(1) = 0.48223;
		u(2) = 0.46120;
		u(3) = 0.45992;
		u(4) = 0.46036;
		u(5) = 0.33556;
		u(6) = 0.20087;
		u(7) = 0.08183;
		u(8) = -0.03039;
		u(9) = -0.07404;
		u(10) = -0.22855;
		u(11) = -0.33050;
		u(12) = -0.40435;
		u(13) = -0.43643;
		u(14) = -0.42901;
		u(15) = -0.41165;
		u(16) = 0.0000;
		return;
	}
	case 7500: {
		u(0) = 1.0;
		u(1) = 0.47244;
		u(2) = 0.47048;
		u(3) = 0.47323;
		u(4) = 0.47167;
		u(5) = 0.34228;
		u(6) = 0.20591;
		u(7) = 0.08342;
		u(8) = -0.03800;
		u(9) = -0.07503;
		u(10) = -0.23176;
		u(11) = -0.32393;
		u(12) = -0.38324;
		u(13) = -0.43025;
		u(14) = -0.43590;
		u(15) = -0.43154;
		u(16) = 0.0000;
		return;
	}
	case 10000: {
		u(0) = 1.0;
		u(1) = 0.47221;
		u(2) = 0.47783;
		u(3) = 0.48070;
		u(4) = 0.47804;
		u(5) = 0.34635;
		u(6) = 0.20673;
		u(7) = 0.08344;
		u(8) = 0.03111;
		u(9) = -0.07540;
		u(10) = -0.23186;
		u(11) = -0.32709;
		u(12) = -0.38000;
		u(13) = -0.41657;
		u(14) = -0.42537;
		u(15) = -0.42735;
		u(16) = 0.0000;
		return;
	}
	default: {
		return;
	}
	}
}

void CompareToGhia::makeReferenceV(numeric_vector& v, size_t reynolds_number) {
	switch (reynolds_number) {
	case 100: {
		v(0) = 0.00000;
		v(1) = -0.05906;
		v(2) = -0.07391;
		v(3) = -0.08864;
		v(4) = -0.10313;
		v(5) = -0.16914;
		v(6) = -0.22445;
		v(7) = -0.24533;
		v(8) = 0.05454;
		v(9) = 0.17527;
		v(10) = 0.17507;
		v(11) = 0.16077;
		v(12) = 0.12317;
		v(13) = 0.10890;
		v(14) = 0.10091;
		v(15) = 0.09233;
		v(16) = 0.00000;
		return;
	}
	case 400: {
		v(0) = 0.0000;
		v(1) = -0.12146;
		v(2) = -0.15663;
		v(3) = -0.19254;
		v(4) = -0.22847;
		v(5) = -0.23827;
		v(6) = -0.44993;
		v(7) = -0.38598;
		v(8) = 0.05186;
		v(9) = 0.30174;
		v(10) = 0.30203;
		v(11) = 0.28124;
		v(12) = 0.22965;
		v(13) = 0.20920;
		v(14) = 0.19713;
		v(15) = 0.18360;
		v(16) = 0.0000;
		return;
	}
	case 1000: {
		v(0) = 0.000;
		v(1) = -0.21388;
		v(2) = -0.27669;
		v(3) = -0.33714;
		v(4) = -0.39188;
		v(5) = -0.51550;
		v(6) = -0.42665;
		v(7) = -0.31966;
		v(8) = 0.02526;
		v(9) = 0.32236;
		v(10) = 0.33075;
		v(11) = 0.37095;
		v(12) = 0.32627;
		v(13) = 0.30353;
		v(14) = 0.29012;
		v(15) = 0.27485;
		v(16) = 0.000;
		return;
	}
	case 3200: {
		v(0) = 0.0000;
		v(1) = -0.39017;
		v(2) = -0.47425;
		v(3) = -0.52357;
		v(4) = -0.54053;
		v(5) = -0.44307;
		v(6) = -0.37401;
		v(7) = -0.31184;
		v(8) = 0.00999;
		v(9) = 0.28188;
		v(10) = 0.29030;
		v(11) = 0.37119;
		v(12) = 0.42768;
		v(13) = 0.41906;
		v(14) = 0.40917;
		v(15) = 0.39560;
		v(16) = 0.00000;
		return;
	}
	case 5000: {
		v(0) = 0.0000;
		v(1) = -0.49774;
		v(2) = -0.55069;
		v(3) = -0.55408;
		v(4) = -0.52876;
		v(5) = -0.41442;
		v(6) = -0.36214;
		v(7) = -0.30018;
		v(8) = 0.00945;
		v(9) = 0.27280;
		v(10) = 0.28066;
		v(11) = 0.35368;
		v(12) = 0.42951;
		v(13) = 0.43648;
		v(14) = 0.43329;
		v(15) = 0.42447;
		v(16) = 0.00000;
		return;
	}
	case 7500: {
		v(0) = 0.0000;
		v(1) = -0.53858;
		v(2) = -0.55216;
		v(3) = -0.52347;
		v(4) = -0.48590;
		v(5) = -0.41050;
		v(6) = -0.36213;
		v(7) = -0.30448;
		v(8) = 0.00824;
		v(9) = 0.27348;
		v(10) = 0.28117;
		v(11) = 0.35060;
		v(12) = 0.41824;
		v(13) = 0.43564;
		v(14) = 0.44030;
		v(15) = 0.43979;
		v(16) = 0.00000;
		return;
	}
	case 10000: {
		v(0) = 0.0000;
		v(1) = -0.54302;
		v(2) = -0.52987;
		v(3) = -0.49099;
		v(4) = -0.45863;
		v(5) = -0.41496;
		v(6) = -0.36737;
		v(7) = -0.30719;
		v(8) = 0.00831;
		v(9) = 0.27224;
		v(10) = 0.28003;
		v(11) = 0.35070;
		v(12) = 0.41487;
		v(13) = 0.43124;
		v(14) = 0.43733;
		v(15) = 0.43983;
		v(16) = 0.00000;
		return;
	}
	default: {
		return;
	}
	}
}

} /* namespace natrium */
