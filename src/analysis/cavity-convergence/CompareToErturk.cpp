/*
 * CompareToErturk.cpp
 *
 *  Created on: 14.03.2016
 *      Author: akraem3m
 */

#include "CompareToErturk.h"

#include <fstream>

#include "deal.II/numerics/fe_field_function.h"
#include "deal.II/base/mpi.h"

namespace natrium {

CompareToErturk::CompareToErturk(const CFDSolver<2>& solver, size_t reynolds,
		string filename) :
		DataProcessor<2>(solver), m_Re(reynolds), m_y(23), m_referenceU(23), m_x(
				23), m_referenceV(23), m_fileName(filename) {
	m_x(0) = 1.0;
	m_x(1) = 0.985;
	m_x(2) = 0.970;
	m_x(3) = 0.955;
	m_x(4) = 0.940;
	m_x(5) = 0.925;
	m_x(6) = 0.910;
	m_x(7) = 0.895;
	m_x(8) = 0.880;
	m_x(9) = 0.865;
	m_x(10) = 0.850;
	m_x(11) = 0.500;
	m_x(12) = 0.150;
	m_x(13) = 0.135;
	m_x(14) = 0.120;
	m_x(15) = 0.105;
	m_x(16) = 0.090;
	m_x(17) = 0.075;
	m_x(18) = 0.060;
	m_x(19) = 0.045;
	m_x(20) = 0.030;
	m_x(21) = 0.015;
	m_x(22) = 0.000;


	m_y(0) = 1.0000;
	m_y(1) = 0.990;
	m_y(2) = 0.980;
	m_y(3) = 0.970;
	m_y(4) = 0.960;
	m_y(5) = 0.950;
	m_y(6) = 0.940;
	m_y(7) = 0.930;
	m_y(8) = 0.920;
	m_y(9) = 0.910;
	m_y(10) = 0.900;
	m_y(11) = 0.500;
	m_y(12) = 0.200;
	m_y(13) = 0.180;
	m_y(14) = 0.160;
	m_y(15) = 0.140;
	m_y(16) = 0.120;
	m_y(17) = 0.100;
	m_y(18) = 0.080;
	m_y(19) = 0.060;
	m_y(20) = 0.040;
	m_y(21) = 0.020;
	m_y(22) = 0.000;
	makeReferenceU(m_referenceU, reynolds);
	makeReferenceV(m_referenceV, reynolds);

}

CompareToErturk::~CompareToErturk() {
	// TODO Auto-generated destructor stub
}

void CompareToErturk::apply() {
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

	// for all test points in Erturk table
	for (size_t i = 0; i < 23; i++) {

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

	m_uError = sqrt(m_uError / 23.);
	m_vError = sqrt(m_vError / 23.);
}

void CompareToErturk::printFinalVelocities() {
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

	// for all test points in Erturk table
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

void CompareToErturk::makeReferenceU(numeric_vector& u, size_t reynolds_number) {
	switch (reynolds_number) {
	case 1000: {
		u(0) = 1.0;
		u(1) = 0.8486;
		u(2) = 0.7065;
		u(3) = 0.5917;
		u(4) = 0.5102;
		u(5) = 0.4582;
		u(6) = 0.4276;
		u(7) = 0.4101;
		u(8) = 0.3993;
		u(9) = 0.3913;
		u(10) = 0.3838;
		u(11) = -0.0620;
		u(12) = -0.3756;
		u(13) = -0.3869;
		u(14) = -0.3854;
		u(15) = -0.3690;
		u(16) = -0.3381;
		u(17) = -0.2960;
		u(18) = -0.2472;
		u(19) = -0.1951;
		u(20) = -0.1392;
		u(21) = -0.0757;
		u(22) = 0.0000;
		return;
	}
	case 2500: {
		u(0) = 1.0;
		u(1) = 0.7704;
		u(2) = 0.5924;
		u(3) = 0.4971;
		u(4) = 0.4607;
		u(5) = 0.4506;
		u(6) = 0.4470;
		u(7) = 0.4424;
		u(8) = 0.4353;
		u(9) = 0.4256;
		u(10) = 0.4141;
		u(11) = -0.0403;
		u(12) = -0.3228;
		u(13) = -0.3439;
		u(14) = -0.3688;
		u(15) = -0.3965;
		u(16) = -0.4200;
		u(17) = -0.4250;
		u(18) = -0.3979;
		u(19) = -0.3372;
		u(20) = -0.2547;
		u(21) = -0.1517;
		u(22) = 0.0000;
		return;
	}
	case 5000: {
		u(0) = 1.0000;
		u(1) = 0.6866;
		u(2) = 0.5159;
		u(3) = 0.4749;
		u(4) = 0.4739;
		u(5) = 0.4738;
		u(6) = 0.4683;
		u(7) = 0.4582;
		u(8) = 0.4452;
		u(9) = 0.4307;
		u(10) = 0.4155;
		u(11) = -0.0319;
		u(12) = -0.3100;
		u(13) = -0.3285;
		u(14) = -0.3467;
		u(15) = -0.3652;
		u(16) = -0.3876;
		u(17) = -0.4168;
		u(18) = -0.4419;
		u(19) = -0.4272;
		u(20) = -0.3480;
		u(21) = -0.2223;
		u(22) = 0.0000;
		return;
	}
	case 7500: {
		u(0) = 1.0000;
		u(1) = 0.6300;
		u(2) = 0.4907;
		u(3) = 0.4817;
		u(4) = 0.4860;
		u(5) = 0.4824;
		u(6) = 0.4723;
		u(7) = 0.4585;
		u(8) = 0.4431;
		u(9) = 0.4275;
		u(10) = 0.4123;
		u(11) = -0.0287;
		u(12) = -0.3038;
		u(13) = -0.3222;
		u(14) = -0.3406;
		u(15) = -0.3587;
		u(16) = -0.3766;
		u(17) = -0.3978;
		u(18) = -0.4284;
		u(19) = -0.4491;
		u(20) = -0.3980;
		u(21) = -0.2633;
		u(22) = 0.0000;
		return;
	}
	case 10000: {
		u(0) = 1.0000;
		u(1) = 0.5891;
		u(2) = 0.4837;
		u(3) = 0.4891;
		u(4) = 0.4917;
		u(5) = 0.4843;
		u(6) = 0.4711;
		u(7) = 0.4556;
		u(8) = 0.4398;
		u(9) = 0.4243;
		u(10) = 0.4095;
		u(11) = -0.0268;
		u(12) = -0.2998;
		u(13) = -0.3179;
		u(14) = -0.3361;
		u(15) = -0.3543;
		u(16) = -0.3721;
		u(17) = -0.3899;
		u(18) = -0.4142;
		u(19) = -0.4469;
		u(20) = -0.4259;
		u(21) = -0.2907;
		u(22) = 0.0000;
		return;
	}
	case 12500: {
		u(0) = 1.0000;
		u(1) = 0.5587;
		u(2) = 0.4833;
		u(3) = 0.4941;
		u(4) = 0.4937;
		u(5) = 0.4833;
		u(6) = 0.4684;
		u(7) = 0.4523;
		u(8) = 0.4366;
		u(9) = 0.4216;
		u(10) = 0.4070;
		u(11) = -0.0256;
		u(12) = -0.2967;
		u(13) = -0.3146;
		u(14) = -0.3326;
		u(15) = -0.3506;
		u(16) = -0.3685;
		u(17) = -0.3859;
		u(18) = -0.4054;
		u(19) = -0.4380;
		u(20) = -0.4407;
		u(21) = -0.3113;
		u(22) = 0.0000;
		return;
	}
	case 15000: {
		u(0) = 1.0000;
		u(1) = 0.5358;
		u(2) = 0.4850;
		u(3) = 0.4969;
		u(4) = 0.4937;
		u(5) = 0.4811;
		u(6) = 0.4653;
		u(7) = 0.4492;
		u(8) = 0.4338;
		u(9) = 0.4190;
		u(10) = 0.4047;
		u(11) = -0.0247;
		u(12) = -0.2942;
		u(13) = -0.3119;
		u(14) = -0.3297;
		u(15) = -0.3474;
		u(16) = -0.3652;
		u(17) = -0.3827;
		u(18) = -0.4001;
		u(19) = -0.4286;
		u(20) = -0.4474;
		u(21) = -0.3278;
		u(22) = 0.0000;
		return;
	}
	case 17500: {
		u(0) = 1.0000;
		u(1) = 0.5183;
		u(2) = 0.4871;
		u(3) = 0.4982;
		u(4) = 0.4925;
		u(5) = 0.4784;
		u(6) = 0.4622;
		u(7) = 0.4463;
		u(8) = 0.4312;
		u(9) = 0.4166;
		u(10) = 0.4024;
		u(11) = -0.0240;
		u(12) = -0.2920;
		u(13) = -0.3096;
		u(14) = -0.3271;
		u(15) = -0.3446;
		u(16) = -0.3622;
		u(17) = -0.3797;
		u(18) = -0.3965;
		u(19) = -0.4206;
		u(20) = -0.4490;
		u(21) = -0.3412;
		u(22) = 0.0000;
		return;
	}
	case 20000: {
		u(0) = 1.0000;
		u(1) = 0.5048;
		u(2) = 0.4889;
		u(3) = 0.4985;
		u(4) = 0.4906;
		u(5) = 0.4754;
		u(6) = 0.4592;
		u(7) = 0.4436;
		u(8) = 0.4287;
		u(9) = 0.4142;
		u(10) = 0.4001;
		u(11) = -0.0234;
		u(12) = -0.2899;
		u(13) = -0.3074;
		u(14) = -0.3248;
		u(15) = -0.3422;
		u(16) = -0.3595;
		u(17) = -0.3769;
		u(18) = -0.3936;
		u(19) = -0.4143;
		u(20) = -0.4475;
		u(21) = -0.3523;
		u(22) = 0.0000;
		return;
	}
	case 21000: {
		u(0) = 1.0000;
		u(1) = 0.5003;
		u(2) = 0.4895;
		u(3) = 0.4983;
		u(4) = 0.4897;
		u(5) = 0.4742;
		u(6) = 0.4580;
		u(7) = 0.4425;
		u(8) = 0.4277;
		u(9) = 0.4132;
		u(10) = 0.3992;
		u(11) = -0.0232;
		u(12) = -0.2892;
		u(13) = -0.3066;
		u(14) = -0.3239;
		u(15) = -0.3412;
		u(16) = -0.3585;
		u(17) = -0.3758;
		u(18) = -0.3925;
		u(19) = -0.4121;
		u(20) = -0.4463;
		u(21) = -0.3562;
		u(22) = 0.0000
		return;
	}
	default: {
		return;
	}
	}
}

void CompareToErturk::makeReferenceV(numeric_vector& v, size_t reynolds_number) {
	switch (reynolds_number) {
	case 1000: {
		v(0) = 0.0000;
		v(1) = -0.0973;
		v(2) = -0.2173;
		v(3) = -0.3400;
		v(4) = -0.4417;
		v(5) = -0.5052;
		v(6) = -0.5263;
		v(7) = -0.5132;
		v(8) = -0.4803;
		v(9) = -0.4407;
		v(1) = -0.4028;
		v(11) = 0.0258;
		v(12) = 0.3756;
		v(13) = 0.3705;
		v(14) = 0.3605;
		v(15) = 0.3460;
		v(16) = 0.3273;
		v(17) = 0.3041;
		v(18) = 0.2746;
		v(19) = 0.2349;
		v(20) = 0.1792;
		v(21) = 0.1019;
		v(22) = 0.0000;
		return;
	}
	case 2500: {
		v(0) = 0.0000;
		v(1) = -0.1675;
		v(2) = -0.3725;
		v(3) = -0.5192;
		v(4) = -0.5603;
		v(5) = -0.5268;
		v(6) = -0.4741;
		v(7) = -0.4321;
		v(8) = -0.4042;
		v(9) = -0.3843;
		v(10) = -0.3671;
		v(11) = 0.0160;
		v(12) = 0.3918;
		v(13) = 0.4078;
		v(14) = 0.4187;
		v(15) = 0.4217;
		v(16) = 0.4142;
		v(17) = 0.3950;
		v(18) = 0.3649;
		v(19) = 0.3238;
		v(20) = 0.2633;
		v(21) = 0.1607;
		v(22) = 0.0000;
		return;
	}
	case 5000: {
		v(0) = 0.0000;
		v(1) = -0.2441;
		v(2) = -0.5019;
		v(3) = -0.5700;
		v(4) = -0.5139;
		v(5) = -0.4595;
		v(6) = -0.4318;
		v(7) = -0.4147;
		v(8) = -0.3982;
		v(9) = -0.3806;
		v(10) = -0.3624;
		v(11) = 0.0117;
		v(12) = 0.3699;
		v(13) = 0.3878;
		v(14) = 0.4070;
		v(15) = 0.4260;
		v(16) = 0.4403;
		v(17) = 0.4426;
		v(18) = 0.4258;
		v(19) = 0.3868;
		v(20) = 0.3263;
		v(21) = 0.2160;
		v(22) = 0.0000;
		return;
	}
	case 7500: {
		v(0) = 0.0000;
		v(1) = -0.2991;
		v(2) = -0.5550;
		v(3) = -0.5434;
		v(4) = -0.4748;
		v(5) = -0.4443;
		v(6) = -0.4283;
		v(7) = -0.4118;
		v(8) = -0.3938;
		v(9) = -0.3755;
		v(10) = -0.3574;
		v(11) = 0.0099;
		v(12) = 0.3616;
		v(13) = 0.3779;
		v(14) = 0.3950;
		v(15) = 0.4137;
		v(16) = 0.4337;
		v(17) = 0.4495;
		v(18) = 0.4494;
		v(19) = 0.4210;
		v(20) = 0.3608;
		v(21) = 0.2509;
		v(22) = 0.0000;
		return;
	}
	case 10000: {
		v(0) = 0.0000;
		v(1) = -0.3419;
		v(2) = -0.5712;
		v(3) = -0.5124;
		v(4) = -0.4592;
		v(5) = -0.4411;
		v(6) = -0.4256;
		v(7) = -0.4078;
		v(8) = -0.3895;
		v(9) = -0.3715;
		v(10) = -0.3538;
		v(11) = 0.0088;
		v(12) = 0.3562;
		v(13) = 0.3722;
		v(14) = 0.3885;
		v(15) = 0.4056;
		v(16) = 0.4247;
		v(17) = 0.4449;
		v(18) = 0.4566;
		v(19) = 0.4409;
		v(20) = 0.3844;
		v(21) = 0.2756;
		v(22) = 0.0000;
		return;
	}
	case 12500: {
		v(0) = 0.0000;
		v(1) = -0.3762;
		v(2) = -0.5694;
		v(3) = -0.4899;
		v(4) = -0.4534;
		v(5) = -0.4388;
		v(6) = -0.4221;
		v(7) = -0.4040;
		v(8) = -0.3859;
		v(9) = -0.3682;
		v(10) = -0.3508;
		v(11) = 0.0080;
		v(12) = 0.3519;
		v(13) = 0.3678;
		v(14) = 0.3840;
		v(15) = 0.4004;
		v(16) = 0.4180;
		v(17) = 0.4383;
		v(18) = 0.4563;
		v(19) = 0.4522;
		v(20) = 0.4018;
		v(21) = 0.2940;
		v(22) = 0.0000;
		return;
	}
	case 15000: {
		v(0) = 0.0000;
		v(1) = -0.4041;
		v(2) = -0.5593;
		v(3) = -0.4754;
		v(4) = -0.4505;
		v(5) = -0.4361;
		v(6) = -0.4186;
		v(7) = -0.4005;
		v(8) = -0.3828;
		v(9) = -0.3654;
		v(10) = -0.3481;
		v(11) = 0.0074;
		v(12) = 0.3483;
		v(13) = 0.3641;
		v(14) = 0.3801;
		v(15) = 0.3964;
		v(16) = 0.4132;
		v(17) = 0.4323;
		v(18) = 0.4529;
		v(19) = 0.4580;
		v(20) = 0.4152;
		v(21) = 0.3083;
		v(22) = 0.0000;
		return;
	}
	case 17500: {
		v(0) = 0.0000;
		v(1) = -0.4269;
		v(2) = -0.5460;
		v(3) = -0.4664;
		v(4) = -0.4482;
		v(5) = -0.4331;
		v(6) = -0.4153;
		v(7) = -0.3975;
		v(8) = -0.3800;
		v(9) = -0.3627;
		v(10) = -0.3457;
		v(11) = 0.0069;
		v(12) = 0.3452;
		v(13) = 0.3608;
		v(14) = 0.3767;
		v(15) = 0.3929;
		v(16) = 0.4093;
		v(17) = 0.4273;
		v(18) = 0.4484;
		v(19) = 0.4602;
		v(20) = 0.4254;
		v(21) = 0.3197;
		v(22) = 0.0000;
		return;
	}
	case 20000: {
		v(0) = 0.0000;
		v(1) = -0.4457;
		v(2) = -0.5321;
		v(3) = -0.4605;
		v(4) = -0.4459;
		v(5) = -0.4300;
		v(6) = -0.4122;
		v(7) = -0.3946;
		v(8) = -0.3774;
		v(9) = -0.3603;
		v(10) = -0.3434;
		v(11) = 0.0065;
		v(12) = 0.3423;
		v(13) = 0.3579;
		v(14) = 0.3736;
		v(15) = 0.3897;
		v(16) = 0.4060;
		v(17) = 0.4232;
		v(18) = 0.4438;
		v(19) = 0.4601;
		v(20) = 0.4332;
		v(21) = 0.3290;
		v(22) = 0.0000;
		return;
	}
	case 21000: {
		v(0) = 0.0000;
		v(1) = -0.4522;
		v(2) = -0.5266;
		v(3) = -0.4588;
		v(4) = -0.4449;
		v(5) = -0.4287;
		v(6) = -0.4110;
		v(7) = -0.3936;
		v(8) = -0.3764;
		v(9) = -0.3593;
		v(10) = -0.3425;
		v(11) = 0.0063;
		v(12) = 0.3413;
		v(13) = 0.3567;
		v(14) = 0.3725;
		v(15) = 0.3885;
		v(16) = 0.4048;
		v(17) = 0.4218;
		v(18) = 0.4420;
		v(19) = 0.4596;
		v(20) = 0.4357;
		v(21) = 0.3323;
		v(22) = 0.0000;
		return;
	}
	default: {
		return;
	}
	}
}

} /* namespace natrium */
