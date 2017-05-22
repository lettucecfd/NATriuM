/*
 * CommandLineParser.cpp
 *
 *  Created on: 02.12.2016
 *      Author: akraem3m
 */

#include "CommandLineParser.h"
#include "boost/type_index.hpp"

namespace natrium {

CommandLineParser::CommandLineParser(int argc, char** argv) :
		m_argc(argc), m_argv(argv), m_description("Allowed options"), m_isImported(
				false) {
	m_description.add_options()("help,h", "produce help message");

	// define standard arguments that can be used in all scripts to change the solver configuration
	setFlag("standard-lbm",
			"sets the order to 2 and cfl to 4*sqrt(2), giving standard lbm on regular grids");
	setArgument<int>("standard-lbm-with",
			"as standard-lbm, but with n-th order calculation of derivatives.");
	setArgument<int>("order-fe", "order of finite elements");
	setArgument<double>("cfl", "CFL number");
	setFlag("output-on", "switches output on");
	setFlag("verbose", "full command line output");
	setArgument<string>("streaming",
			"advection scheme (sedg (discontinuous Galerkin) or sl (semi-Lagrange))");
	setArgument<string>("collision",
			"collision scheme (bgk, bgkt (transformed), inc (incompressible), pp (pseudo-potential)"
					", ss (steady state), mrt, kbcs (KBC standard), kbcc (KBC central), reg (regularized by Latt)"
					", am4 (multistep Adams-Moulton 4), bdf2 (multistep BDF2), est (entropic stablized)");
	setArgument<string>("integrator",
			"time integrator in sedg streaming ["
					"(1) rk4 (classical Runge-Kutta, built-in), (2) theta (theta method), (3) exp (exponential),  (4) ee (explicit Euler), "
					"(5) rk3 (explicit 3rd order RK), (6) rk4_2 (classical RK, deal.ii), (7) ie (implicit Euler), (8) midpoint (implicit midpoint), "
					"(9) cn (crank-nicolson), (10) sdirk (diagonally-implicit RK with two stages), (11) heun-euler (Heun-Euler), (12) bogacki (Bogacki-Shampine), "
					"(13) dopri (Dorman-Price), (14) fehlberg (Fehlberg), (15) cash-karp (Cash-Karp)]");
	setArgument<int>("integrator-id",
			"select integrator by id, see help message for --integrator");
	setArgument<double>("scaling", "scaling of discrete particle velocities");
	setArgument<string>("regularization",
			"regularization of high-order moments [no (no regularization), pem (pseudo-entropy maximization), "
					"zero (zero high-order moments), em (entropy maximization), pemwe(pseudo-entropy maximization with flexible moment e)]");
	setArgument<string>("support-points",
			"support points for semi-Lagrangian streaming [gll (Gauss-Lobatto-Legendre)"
					", glc (Gauss-Lobatto-Chebyshev), gc (Gauss-Chebyshev), equi (equidistant)]");
	setArgument<string>("output-dir", "output directory");
	setArgument<int>("output-sol", "output solution interval (#iterations)");
	setArgument<int>("output-chk", "output checkpoint interval (#iterations)");
	setArgument<int>("output-tab", "output table interval (#iterations)");
	setArgument<string>("stencil",
			"stencil that defines the discrete particle velocities [d2q9, d3q19, d3q15, d3q27]");
	setArgument<double>("tmax", "simulation end time");
	setArgument<string>("init",
			"initialization scheme [equi (equilibrium), iter (iterative)]");
	setArgument<int>("init-niter",
			"iterative initialization scheme: max. number of iterations");
	setArgument<double>("init-res",
			"iterative initialization scheme: residual");

}

void CommandLineParser::applyToSolverConfiguration(SolverConfiguration& cfg) {

	// standard lbm on regular grids
	if (hasArgument("standard-lbm")) {
		if (hasArgument("cfl")) {
			throw CommandLineParserException(
					"Conflicting arguments: standard-lbm and cfl");
		}
		if (hasArgument("order-fe") and (getArgument<int>("order-fe") != 2)) {
			throw CommandLineParserException(
					"Conflicting arguments: standard-lbm and order-fe != 2");
		}
		if (hasArgument("streaming")
				and (getArgument<string>("streaming") != "sl")) {
			throw CommandLineParserException(
					"Conflicting arguments: standard-lbm and streaming != sl");
		}
		cfg.setCFL(sqrt(2) * 2);
		cfg.setSedgOrderOfFiniteElement(2);
		cfg.setAdvectionScheme(SEMI_LAGRANGIAN);
		LOG(BASIC) << "Scheme set to standard LBM via command line." << endl;
	}

	// standard lbm on regular grids with higher-order differentiation
	if (hasArgument("standard-lbm-with")) {
		if (hasArgument("cfl")) {
			throw CommandLineParserException(
					"Conflicting arguments: standard-lbm and cfl");
		}
		if (hasArgument("order-fe") and (getArgument<int>("order-fe") != 2)) {
			throw CommandLineParserException(
					"Conflicting arguments: standard-lbm and order-fe != 2");
		}
		if (hasArgument("streaming")
				and (getArgument<string>("streaming") != "sl")) {
			throw CommandLineParserException(
					"Conflicting arguments: standard-lbm and streaming != sl");
		}
		if (hasArgument("standard-lbm")) {
			throw CommandLineParserException(
					"Conflicting arguments: standard-lbm and standard-lbm-with");
		}
		int order = getArgument<int>("standard-lbm-with");
		assert(order > 0);
		cfg.setCFL(sqrt(2) * order * order / order);
		cfg.setSedgOrderOfFiniteElement(order);
		cfg.setSupportPoints(EQUIDISTANT_POINTS);
		cfg.setAdvectionScheme(SEMI_LAGRANGIAN);
		LOG(BASIC) << "Scheme set to standard LBM with N="
				<< cfg.getSedgOrderOfFiniteElement() << " via command line."
				<< endl;
	}

	// finite element order
	if (hasArgument("order-fe")) {
		int order_fe = getArgument<int>("order-fe");
		cfg.setSedgOrderOfFiniteElement(order_fe);
		LOG(BASIC) << "Order set to " << order_fe << " via command line."
				<< endl;
	}

	// cfl number
	if (hasArgument("cfl")) {
		double cfl = getArgument<double>("cfl");
		cfg.setCFL(cfl);
		LOG(BASIC) << "CFL set to " << cfl << " via command line." << endl;
	}

	// switch output on
	if (hasArgument("output-on")) {
		cfg.setSwitchOutputOff(false);
		LOG(BASIC) << "Output switched on via command line." << endl;
	}

	// command line verbosity
	if (hasArgument("verbose")) {
		cfg.setCommandLineVerbosity(static_cast<int>(ALL));
		cfg.setSwitchOutputOff(false);
		LOG(BASIC) << "Max verbosity specified via command line." << endl;
	}

	// advection scheme
	if (hasArgument("streaming")) {
		string streaming = getArgument<string>("streaming");
		if (streaming == "sedg") {
			cfg.setAdvectionScheme(SEDG);
		} else if (streaming == "sl") {
			cfg.setAdvectionScheme(SEMI_LAGRANGIAN);
		} else {
			throw CommandLineParserException(
					"--streaming had illegal value (allowed values: sl and sedg)");
		}
		LOG(BASIC) << "Advection scheme set to " << streaming
				<< " via command line." << endl;
	}

	// collision scheme
	if (hasArgument("collision")) {
		string collision = getArgument<string>("collision");
		if (collision == "bgk")
			cfg.setCollisionScheme(BGK_STANDARD);
		else if (collision == "bgkt")
			cfg.setCollisionScheme(BGK_STANDARD_TRANSFORMED);
		else if (collision == "inc")
			cfg.setCollisionScheme(BGK_INCOMPRESSIBLE);
		else if (collision == "pp")
			cfg.setCollisionScheme(BGK_MULTIPHASE);
		else if (collision == "ss")
			cfg.setCollisionScheme(BGK_STEADY_STATE);
		else if (collision == "mrt")
			cfg.setCollisionScheme(MRT_STANDARD);
		else if (collision == "kbcs")
			cfg.setCollisionScheme(KBC_STANDARD);
		else if (collision == "kbcc")
			cfg.setCollisionScheme(KBC_CENTRAL);
		else if (collision == "am4")
			cfg.setCollisionScheme(BGK_MULTI_AM4);
		else if (collision == "bdf2")
			cfg.setCollisionScheme(BGK_MULTI_BDF2);
		else if (collision == "reg")
			cfg.setCollisionScheme(BGK_REGULARIZED);
		else if (collision == "est")
			cfg.setCollisionScheme(ENTROPIC_STABILIZED);
		else {
			throw CommandLineParserException(
					"--collision had illegal value (see --help)");
		}
		LOG(BASIC) << "Collision set to " << collision << " via command line."
				<< endl;
	}

	// time integrator
	if (hasArgument("integrator")) {
		if (hasArgument("integrator-id")) {
			throw CommandLineParserException(
					"Conflicting arguments: integrator and integrator-id");
		}
		string integrator = getArgument<string>("integrator");
		if (integrator == "rk4") {
			cfg.setTimeIntegrator(RUNGE_KUTTA_5STAGE);
			cfg.setDealIntegrator(NONE);
		} else if (integrator == "theta") {
			cfg.setTimeIntegrator(THETA_METHOD);
			cfg.setDealIntegrator(NONE);
		} else if (integrator == "exp") {
			cfg.setTimeIntegrator(EXPONENTIAL);
			cfg.setDealIntegrator(NONE);
		} else if (integrator == "ee") {
			cfg.setTimeIntegrator(OTHER);
			cfg.setDealIntegrator(FORWARD_EULER);
		} else if (integrator == "rk3") {
			cfg.setTimeIntegrator(OTHER);
			cfg.setDealIntegrator(RK_THIRD_ORDER);
		} else if (integrator == "rk4_2") {
			cfg.setTimeIntegrator(OTHER);
			cfg.setDealIntegrator(RK_CLASSIC_FOURTH_ORDER);
		} else if (integrator == "ie") {
			cfg.setTimeIntegrator(OTHER);
			cfg.setDealIntegrator(BACKWARD_EULER);
		} else if (integrator == "midpoint") {
			cfg.setTimeIntegrator(OTHER);
			cfg.setDealIntegrator(IMPLICIT_MIDPOINT);
		} else if (integrator == "cn") {
			cfg.setTimeIntegrator(OTHER);
			cfg.setDealIntegrator(CRANK_NICOLSON);
		} else if (integrator == "sdirk") {
			cfg.setTimeIntegrator(OTHER);
			cfg.setDealIntegrator(SDIRK_TWO_STAGES);
		} else if (integrator == "heun-euler") {
			cfg.setTimeIntegrator(OTHER);
			cfg.setDealIntegrator(HEUN_EULER);
		} else if (integrator == "bogacki") {
			cfg.setTimeIntegrator(OTHER);
			cfg.setDealIntegrator(BOGACKI_SHAMPINE);
		} else if (integrator == "dopri") {
			cfg.setTimeIntegrator(OTHER);
			cfg.setDealIntegrator(DOPRI);
		} else if (integrator == "fehlberg") {
			cfg.setTimeIntegrator(OTHER);
			cfg.setDealIntegrator(FEHLBERG);
		} else if (integrator == "cash-karp") {
			cfg.setTimeIntegrator(OTHER);
			cfg.setDealIntegrator(CASH_KARP);
		} else {
			throw CommandLineParserException(
					"--integrator had illegal value (see --help)");
		}
		LOG(BASIC) << "Time integrator set to " << integrator
				<< " via command line" << endl;
	}

	// integrator-id
	if (hasArgument("integrator-id")) {
		if (hasArgument("integrator")) {
			throw CommandLineParserException(
					"Conflicting arguments: integrator and integrator-id");
		}
		int id = getArgument<int>("integrator-id");
		TimeIntegratorName t;
		DealIntegratorName d;
		string s;
		CFDSolverUtilities::get_integrator_by_id(id, t, d, s);
		cfg.setTimeIntegrator(t);
		cfg.setDealIntegrator(d);
		LOG(BASIC) << "Time integrator set to " << s << " via command line"
				<< endl;
	}

	// stencil scaling
	if (hasArgument("scaling")) {
		double scaling = getArgument<double>("scaling");
		cfg.setStencilScaling(scaling);
		LOG(BASIC) << "Scaling set to " << scaling << " via command line"
				<< endl;
	}

	// regularization scheme
	if (hasArgument("regularization")) {
		string reg = getArgument<string>("regularization");
		if (reg == "no") {
			cfg.setRegularizationScheme(NO_REGULARIZATION);
		} else if (reg == "pem") {
			cfg.setRegularizationScheme(PSEUDO_ENTROPY_MAXIMIZATION);
		} else if (reg == "zero") {
			cfg.setRegularizationScheme(ZERO_HIGH_ORDER_MOMENTS);
		} else if (reg == "em") {
			cfg.setRegularizationScheme(ENTROPY_MAXIMIZATION);
		} else if (reg == "pemwe") {
			cfg.setRegularizationScheme(PSEUDO_ENTROPY_MAXIMIZATION_WITH_E);
		} else {
			std::stringstream msg;
			msg << "Regularization scheme " << reg << " is illegal." << endl;
			msg << "See --help for allowed options." << endl;
			throw CommandLineParserException(msg.str());
		}
		LOG(BASIC) << "Regularization set to " << reg << " via command line"
				<< endl;
	}

	// support points
	if (hasArgument("support-points")) {
		string sup_p = getArgument<string>("support-points");
		if (sup_p == "gll") {
			cfg.setSupportPoints(GAUSS_LOBATTO_POINTS);
		} else if (sup_p == "glc") {
			cfg.setSupportPoints(GAUSS_LOBATTO_CHEBYSHEV_POINTS);
		} else if (sup_p == "gc") {
			cfg.setSupportPoints(GAUSS_CHEBYSHEV_POINTS);
		} else if (sup_p == "equi") {
			cfg.setSupportPoints(EQUIDISTANT_POINTS);
		} else {
			std::stringstream msg;
			msg << "--support-points=" << sup_p << " is illegal." << endl;
			msg << "See --help for allowed options." << endl;
			throw CommandLineParserException(msg.str());
		}
		LOG(BASIC) << "Support points set to " << sup_p << " via command line"
				<< endl;
	}

	// output directory
	if (hasArgument("output-dir")) {
		string dir = getArgument<string>("output-dir");
		cfg.setOutputDirectory(dir);
		cfg.setSwitchOutputOff(false);
		LOG(BASIC) << "Output directory set to " << dir << " via command line"
				<< endl;
	}

	// output solution interval
	if (hasArgument("output-sol")) {
		int sol = getArgument<int>("output-sol");
		cfg.setOutputSolutionInterval(sol);
		cfg.setSwitchOutputOff(false);
		LOG(BASIC) << "Output solution interval set to " << sol
				<< " via command line" << endl;
	}

	// output checkpoint interval
	if (hasArgument("output-chk")) {
		int erval = getArgument<int>("output-chk");
		cfg.setOutputCheckpointInterval(erval);
		cfg.setSwitchOutputOff(false);
		LOG(BASIC) << "Output checkpoint interval set to " << erval
				<< " via command line" << endl;
	}

	// output solution interval
	if (hasArgument("output-tab")) {
		int erval = getArgument<int>("output-tab");
		cfg.setOutputTableInterval(erval);
		cfg.setSwitchOutputOff(false);
		LOG(BASIC) << "Output table interval set to " << erval
				<< " via command line" << endl;
	}

	// stencil
	if (hasArgument("stencil")) {
		const string sten = getArgument<string>("stencil");
		if (sten == "d2q9") {
			cfg.setStencil(Stencil_D2Q9);
		} else if (sten == "d3q19") {
			cfg.setStencil(Stencil_D3Q19);
		} else if (sten == "d3q15") {
			cfg.setStencil(Stencil_D3Q15);
		} else if (sten == "d3q27") {
			cfg.setStencil(Stencil_D3Q27);
		} else {
			std::stringstream msg;
			msg << "--stencil=" << sten << " is illegal." << endl;
			msg << "See --help for allowed options." << endl;
			throw CommandLineParserException(msg.str());
		}
		LOG(BASIC) << "Stencil set to " << sten << " via command line" << endl;
	}

	// simulation end time
	if (hasArgument("tmax")) {
		int tmax = getArgument<double>("tmax");
		cfg.setSimulationEndTime(tmax);
		LOG(BASIC) << "Output simulation end time set to " << tmax
				<< " via command line" << endl;
	}

	// iterative initialization
	if (hasArgument("init")) {
		const string init = getArgument<string>("init");
		if (init == "equi") {
			cfg.setInitializationScheme(EQUILIBRIUM);
		} else if (init == "iter") {
			cfg.setInitializationScheme(ITERATIVE);
		} else {
			std::stringstream msg;
			msg << "--init=" << init << " is illegal." << endl;
			msg << "See --help for allowed options." << endl;
			throw CommandLineParserException(msg.str());
		}
		LOG(BASIC) << "Initialization set to " << init << " via command line" << endl;
	}

	// iterative initialization: residual
	if (hasArgument("init-res")) {
		const double res = getArgument<double>("init-res");
		cfg.setIterativeInitializationResidual(res);
		LOG(BASIC) << "Residual for iterative initialization set to " << res << " via command line" << endl;
	}

	// iterative initialization: max n. iterations
	if (hasArgument("init-niter")) {
		const int niter = getArgument<int>("init-niter");
		cfg.setIterativeInitializationNumberOfIterations(niter);
		LOG(BASIC) << "Max. number of iterative initialization steps set to " << niter << " via command line" << endl;
	}

}

template<class type>
void CommandLineParser::setPositionalArgument(std::string name,
		std::string description) {
	assert_clean();
	try {
		m_description.add_options()(name.c_str(), po::value<type>(),
				description.c_str());
		m_positionalOptions.add(name.c_str(), 1);
	} catch (po::error& e) {
		std::stringstream msg;
		msg << "A problem occured when defining the positional argument ";
		msg << name;
		msg << "Make sure that the name of the argument is unique." << endl;
		msg << endl;
		msg << "Original error message:" << e.what();
		msg << "The following arguments were set previously:" << endl;
		msg << m_description << endl;
		throw CommandLineParserException(msg.str());
	}
}

void CommandLineParser::makePositional(std::string name) {
	assert_clean();
	/*if (not hasArgument(name)){
	 std::stringstream msg;
	 msg << "A problem occured when making ";
	 msg << name;
	 msg << " a positional argument." << endl;
	 msg << endl;
	 msg << "The argument '" << name << "' was not specified before the call to makePositional()." << endl;
	 msg << m_description << endl;
	 throw CommandLineParserException(msg.str());

	 }*/
	try {
		m_positionalOptions.add(name.c_str(), 1);
	} catch (po::error& e) {
		std::stringstream msg;
		msg << "A problem occured when defining the positional argument ";
		msg << name;
		msg
				<< "Make sure that the name of the argument was defined (but not as a positional argument) before."
				<< endl;
		msg << endl;
		msg << "Original error message:" << e.what();
		msg << "Help message:" << endl;
		msg << *this << endl;
		throw CommandLineParserException(msg.str());
	}
}

//explicit instantiation
template
void CommandLineParser::setPositionalArgument<int>(std::string name,
		std::string description);
template
void CommandLineParser::setPositionalArgument<double>(std::string name,
		std::string description);
template
void CommandLineParser::setPositionalArgument<string>(std::string name,
		std::string description);

template<class type>
void CommandLineParser::setArgument(std::string name, std::string description,
		type default_value) {
	assert_clean();
	try {
		type opt;
		m_description.add_options()(name.c_str(),
				po::value<type>(&opt)->default_value(default_value),
				description.c_str());
	} catch (po::error& e) {
		std::stringstream msg;
		msg << "A problem occured when defining the argument ";
		msg << name;
		msg << "Make sure that the name of the argument is unique." << endl;
		msg << endl;
		msg << "Original error message:" << e.what();
		msg << "The following arguments were set previously:" << endl;
		msg << m_description << endl;
		throw CommandLineParserException(msg.str());
	}
}

template void CommandLineParser::setArgument<int>(std::string name,
		std::string description, int default_value);
template void CommandLineParser::setArgument<double>(std::string name,
		std::string description, double default_value);
template void CommandLineParser::setArgument<string>(std::string name,
		std::string description, string default_value);

template<class type>
void CommandLineParser::setArgument(std::string name, std::string description) {
	assert_clean();
	try {
		m_description.add_options()(name.c_str(), po::value<type>(),
				description.c_str());
	} catch (po::error& e) {
		std::stringstream msg;
		msg << "A problem occured when defining the argument ";
		msg << name;
		msg << "Make sure that the name of the argument is unique." << endl;
		msg << endl;
		msg << "Original error message:" << e.what();
		msg << "The following arguments were set previously:" << endl;
		msg << m_description << endl;
		throw CommandLineParserException(msg.str());
	}
}
template void CommandLineParser::setArgument<int>(std::string name,
		std::string description);
template void CommandLineParser::setArgument<double>(std::string name,
		std::string description);
template void CommandLineParser::setArgument<string>(std::string name,
		std::string description);

template<class type>
type CommandLineParser::getArgument(std::string name) {
	assert_imported();
	if (m_varMap.count(name)) {
		try {
			return m_varMap[name].as<type>();
		} catch (po::error& e) {
			std::stringstream msg;
			msg << "A problem occured when reading the argument ";
			msg << name;
			msg << " (of type double.)"
					<< "Make sure that the argument is properly defined.";
			msg << endl;
			msg << "Original error message:" << e.what();
			throw CommandLineParserException(msg.str());
		} catch (boost::exception_detail::clone_impl<
				boost::exception_detail::error_info_injector<boost::bad_any_cast> >& e) {
			std::stringstream msg;
			msg << "A problem occured when reading the argument ";
			msg << name;
			msg << " with getArgument<";
			msg << boost::typeindex::type_id<type>().pretty_name();
			msg << ">. The argument was defined with a different type." << endl;
			msg << "Original error message:" << e.what();
			throw CommandLineParserException(msg.str());
		}
	} else {
		std::stringstream msg;
		msg << "Did not find command line argument " << name;
		msg << ". Pass --help to see the help message." << endl;
		throw CommandLineParserException(msg.str());
	}
}
template double CommandLineParser::getArgument<double>(std::string name);
template int CommandLineParser::getArgument<int>(std::string name);
template string CommandLineParser::getArgument<string>(std::string name);

void CommandLineParser::setFlag(std::string name, std::string description) {
	assert_clean();
	try {
		m_description.add_options()(name.c_str(), description.c_str());
	} catch (po::error& e) {
		std::stringstream msg;
		msg << "A problem occured when defining the flag ";
		msg << name;
		msg << " (of type flag.)"
				<< "Make sure that the name of the argument is unique." << endl;
		msg << endl;
		msg << "Original error message:" << e.what();
		throw CommandLineParserException(msg.str());
	}
}

void CommandLineParser::importOptions() {

	m_isImported = true;

	// positional arguments
	po::store(
			po::command_line_parser(m_argc, m_argv).options(m_description).positional(
					m_positionalOptions).run(), m_varMap);
	po::notify(m_varMap);

	// print out help message and stop program
	if (m_varMap.count("help")) {
		pout << *this << endl;
		throw HelpMessageStop();
	}
}

std::ostream& operator<<(std::ostream& os, const CommandLineParser& obj) {
	os << endl;
	os << "Help message for NATriuM executable " << obj.m_name << endl;
	os << "------------------------------------------------------" << endl;
	os << "Usage: \t" << obj.m_name << "  ";
	for (size_t i = 0; i < obj.m_positionalOptions.max_total_count(); i++) {
		os << "<" << obj.m_positionalOptions.name_for_position(i) << "> ";
	}
	os << endl;
	os << endl;
	os << "Documentation:" << endl;
	os << obj.m_documentation << endl;
	os << "NATriuM is a library for off-lattice Boltzmann simulations." << endl;
	os << endl;
	os << obj.m_description << endl;
	return os;
}

} /* namespace natrium */

