/**
 * @file SEDGMinLee_test.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "advection/SEDGMinLee.h"

#include <fstream>
#include <map>
#include "complex.h"
#include <cstring>
#include <string>
#include <dirent.h>
#include <sys/stat.h>
#include <stdlib.h>

#include "boost/test/unit_test.hpp"

#include "deal.II/numerics/data_out.h"
#include "deal.II/fe/fe_dgq.h"
#include "deal.II/fe/fe_update_flags.h"
#include "deal.II/dofs/dof_tools.h"
#include "deal.II/base/quadrature_lib.h"
#include "deal.II/lac/lapack_full_matrix.h"

#include "boltzmannmodels/D2Q9IncompressibleModel.h"

#include "timeintegration/RungeKutta5LowStorage.h"

#include "PeriodicTestDomain2D.h"

#include "utilities/BasicNames.h"

///////////////////////////////
///////////////////////////////
// FLAGS TO SWITCH OUTPUT ON //
 #define CREATE_DATA_FILES // This takes a lot of time (~1h)
 #define PRINT_SYSTEM_MATRIX
 #define EULER_OUT
 #define RK5_OUT
///////////////////////////////
///////////////////////////////

namespace natrium {

BOOST_AUTO_TEST_SUITE(SEDGMinLee_test)

BOOST_AUTO_TEST_CASE(SEDGMinLee_Construction_test) {
	cout << "SEDGMinLee_Construction_test..." << endl;

	size_t fe_order = 2;
	size_t refinementLevel = 2;
	PeriodicTestDomain2D periodic(refinementLevel);
	BOOST_CHECK_NO_THROW(
			SEDGMinLee<2> streaming(periodic.getTriangulation(), periodic.getBoundaries(), fe_order, make_shared<D2Q9IncompressibleModel>()));

	cout << "done." << endl;
} /* SEDGMinLee_Construction_test */

BOOST_AUTO_TEST_CASE(SEDGMinLee_systemMatrix_test) {
	cout << "SEDGMinLee_systemMatrix_test..." << endl;

	bool useLaxFlux = true;
	do { // useCentralFlux = true/false
#ifdef CREATE_DATA_FILES
	// Make results dir
		std::stringstream dirname;
		dirname << getenv("NATRIUM_HOME") << "/maxeigenvalue_advection";
	if (mkdir(dirname.str().c_str(), 0777) == -1) {
		if (errno != EEXIST) {
			cerr << "Fehler in mkdir: " << strerror(errno) << endl;
		}
	}

	// Make results dir 2
	std::stringstream dirname2;
	dirname2 << getenv("NATRIUM_HOME") << "/eigenvalues_advection";
	if (mkdir(dirname2.str().c_str(), 0777) == -1) {
		if (errno != EEXIST) {
			cerr << "Fehler in mkdir: " << strerror(errno) << endl;
		}
	}

	for (size_t fe_order = 3; fe_order <= 9; fe_order+=2) {
		// create files for max eigenvalue diagram
		std::stringstream maxEigenvalueFilename;
		maxEigenvalueFilename << dirname;
		if (!useLaxFlux) {
			maxEigenvalueFilename << "/maxeigenvalue_cflux_";
		} else {
			maxEigenvalueFilename << "/maxeigenvalue_laxflux_";
		}
		maxEigenvalueFilename << "feorder_" << fe_order << ".dat";
		cout << maxEigenvalueFilename.str().c_str() << "..." << endl;
		std::ofstream normOut(maxEigenvalueFilename.str().c_str());
#else
		for (size_t fe_order = 2; fe_order < 3; fe_order++) {
#endif
			for (size_t refinementLevel = 1; refinementLevel < 3;
					refinementLevel++) {
				PeriodicTestDomain2D periodic(refinementLevel);
				SEDGMinLee<2> streaming(periodic.getTriangulation(),
						periodic.getBoundaries(), fe_order,
						make_shared<D2Q9IncompressibleModel>(), "",
						!useLaxFlux);
				const distributed_sparse_block_matrix& matrices =
						streaming.getSystemMatrix();
				BOOST_CHECK(matrices.n_block_cols() == D2Q9IncompressibleModel::Q-1);
				BOOST_CHECK(matrices.n_block_rows() == D2Q9IncompressibleModel::Q-1);
#ifdef PRINT_SYSTEM_MATRIX
				for (size_t i = 0; i < 2; i++) { //D2Q9IncompressibleModel::Q; i++){
					cout << "Matrix " << i << ": " << endl;
					matrices.block(i,i).print_formatted(cout);
				}
#endif
				//BOOST_CHECK()
				//singular value decomposition
				dealii::LAPACKFullMatrix<double> systemMatrix(
						matrices.block(0,0).n());
				for (size_t i = 0; i < D2Q9IncompressibleModel::Q-1; i++) {
#ifdef CREATE_DATA_FILES
					// create files for spectrum plots
					std::stringstream filename;
					filename << dirname2.str();
					if (!useLaxFlux) {
						filename << "/cflux_";
					} else {
						filename << "/laxflux_";
					}
					filename << "refine_" << refinementLevel << "_order_"
					<< fe_order << "_matrix_" << i << "_spectrum.dat";
					cout << filename.str().c_str() << "..." << endl;
					std::ofstream spectrumOut(filename.str().c_str());
#endif
					// compute eigenvalues
					systemMatrix = 0;
					double deltaX = 1.
							/ (pow(2, refinementLevel) * (fe_order - 1));
					double deltaT = deltaX;
					distributed_sparse_matrix sparseSytemMatrix(
							matrices.block(i,i).get_sparsity_pattern());
					sparseSytemMatrix.copy_from(matrices.block(i,i));
					sparseSytemMatrix *= deltaT;
					systemMatrix.copy_from(sparseSytemMatrix);
					systemMatrix.compute_eigenvalues();
					// check eigenvalues (the real part must vanish for central flux, the absolute value must be in the same order of magnitude as delta x)
					double max_eigenvalue = 0;
					for (size_t j = 0; j < systemMatrix.n_cols(); j++) {
						std::
						complex<double> eigenValue;
						eigenValue = systemMatrix.eigenvalue(j);
						BOOST_CHECK(eigenValue.real() < 0.1);
						//BOOST_CHECK(eigenValue.imag() < 1.0);
						if (!useLaxFlux) {
							BOOST_CHECK(fabs(eigenValue.real()) < 1e-10);
						}
						if (abs(eigenValue) > max_eigenvalue) {
							max_eigenvalue = abs(eigenValue);
						}
					}
					if (i != 0) {
						BOOST_CHECK(max_eigenvalue > 0.5 * deltaX);
						//BOOST_CHECK(max_eigenvalue < 2.5 * deltaX);
					}
#ifdef CREATE_DATA_FILES
					// print out eigenvalues and spectralnorm
					for (size_t j = 0; j < systemMatrix.n_cols(); j++) {
						std::
						complex<double> eigenValue;
						eigenValue = systemMatrix.eigenvalue(j);
						spectrumOut << eigenValue.real() << " "
						<< eigenValue.imag() << endl;
					}
					if (i == 1) {
						normOut << refinementLevel << " " << max_eigenvalue
						<< endl;
					}
#endif
				}
			}
		}
		useLaxFlux = !useLaxFlux;
	} while (!useLaxFlux);

	cout << "done." << endl;
}
/* SEDGMinLee_systemMatrix_test */

BOOST_AUTO_TEST_CASE(SEDGMinLee_steadyStreaming_test) {
	cout << "SEDGMinLee_steadyStreaming_test..." << endl;

	size_t refinementLevel = 3;
	size_t fe_order = 3;
	PeriodicTestDomain2D periodic(refinementLevel);

	SEDGMinLee<2> streaming(periodic.getTriangulation(),
			periodic.getBoundaries(), fe_order,
			make_shared<D2Q9IncompressibleModel>());
	const distributed_sparse_block_matrix& matrices =
			streaming.getSystemMatrix();
	// choose time step and number of time steps so dx = dt and the bump passes the domain one time
	const double timeStep = 1. / (pow(2, refinementLevel) * (fe_order - 1));
	const size_t numberOfTimeSteps = size_t(1. / timeStep);

	// Initialize all particle distribution functions with 1
	distributed_vector f(streaming.getNumberOfDoFs());
	for (size_t i = 0; i < f.size(); i++) {
		f(i) = 1.0;
	}

	// CHECK MASS CONSERVATION AND CONSERVATION OF STEADY STATE
	double initialMass = f.size();
	// Explicit euler for streaming in direction (1,0)
	distributed_sparse_matrix advectionMatrix;
	advectionMatrix.reinit(streaming.getSparsityPattern(2));
	advectionMatrix.copy_from(matrices.block(2,2));
	advectionMatrix *= timeStep;
	distributed_vector f_tmp(streaming.getNumberOfDoFs());
	for (size_t i = 0; i < numberOfTimeSteps; i++) {
		//stream
		f_tmp = f;
		advectionMatrix.vmult_add(f, f_tmp);
		double mass = 0.0;
		for (size_t j = 0; j < f.size(); j++) {
			BOOST_CHECK(fabs(f(j) - 1.0) < 1e-10);
			mass += f(j);
		}
		BOOST_CHECK(fabs(mass - initialMass) / initialMass < 1e-5);

	}

	cout << "done." << endl;
} /* SEDGMinLee_steadyStreaming_test */

BOOST_AUTO_TEST_CASE(SEDGMinLee_streaming_test) {
	cout << "SEDGMinLee_streaming_test..." << endl;

	size_t refinementLevel = 3;
	size_t fe_order = 3;
	PeriodicTestDomain2D periodic(refinementLevel);

	SEDGMinLee<2> streaming(periodic.getTriangulation(),
			periodic.getBoundaries(), fe_order,
			make_shared<D2Q9IncompressibleModel>());
	const distributed_sparse_block_matrix& matrices =
			streaming.getSystemMatrix();
	// choose time step and number of time steps so dx = dt and the bump passes the domain one time
	const double timeStep = 0.1 / (pow(2, refinementLevel) * (fe_order - 1));
	const size_t numberOfTimeSteps = size_t(1. / timeStep);

	// Initialize all particle distribution functions with 1, one corner element with 2
	distributed_vector f(streaming.getNumberOfDoFs());
	double initialMass = 0.0;
	for (size_t i = 0; i < f.size(); i++) {
		f(i) = 1.0;
		initialMass += 1;
	}
	// create smooth initial conditions
	vector<dealii::Point<2> > supportPoints(
			streaming.getDoFHandler()->n_dofs());
	dealii::DoFTools::map_dofs_to_support_points(streaming.getMapping(),
			*streaming.getDoFHandler(), supportPoints);
	dealii::Point<2> midPoint(0.25, 0.25);
	const double PI = 3.14;
	for (size_t i = 0; i < streaming.getDoFHandler()->n_dofs(); i++) {
		double distance = supportPoints.at(i).distance(midPoint);
		if (distance <= 0.25) {
			f(i) = 1 + cos(PI / 0.5 * distance);
			initialMass += cos(PI / 0.5 * distance);
		}
	}
	// CHECK MASS CONSERVATION
	// Explicit euler for streaming in direction (1,0)
	distributed_sparse_matrix advectionMatrix;
	advectionMatrix.reinit(streaming.getSparsityPattern(2));
	advectionMatrix.copy_from(matrices.block(2,2));
	advectionMatrix *= timeStep;
	distributed_vector f_tmp(streaming.getNumberOfDoFs());

	// Make results dir
#ifdef EULER_OUT
	std::stringstream dirname;
	dirname << getenv("NATRIUM_HOME") << "/expliciteuler_periodicstreaming";
	if (mkdir(dirname.str().c_str(), 0777) == -1) {
		if (errno != EEXIST) {
			cerr << "Fehler in mkdir: " << strerror(errno) << endl;
		}
	}
#endif
	for (size_t i = 0; i < numberOfTimeSteps; i++) {
#ifdef EULER_OUT
		// output
		std::stringstream str;
		str << dirname.str() << "/t_" << i << ".dat";
		std::string filename = str.str();
		std::ofstream gnuplot_output (filename.c_str());
		dealii::DataOut<2> data_out;
		data_out.attach_dof_handler (*streaming.getDoFHandler());
		data_out.add_data_vector (f, "f");
		data_out.build_patches ();
		data_out.write_gnuplot(gnuplot_output);
#endif
		//stream
		f_tmp = f;
		advectionMatrix.vmult_add(f, f_tmp);

		double mass = 0.0;
		for (size_t j = 0; j < f.size(); j++) {
			mass += f(j);
		}

		BOOST_CHECK(fabs(mass - initialMass) / initialMass < 5e-2);
		// NOTE: mass is not conserved exactly for explicit euler (but up to 5 percent)

	}
	double mass = 0.0;
	for (size_t j = 0; j < f.size(); j++) {
		mass += f(j);
	}
	BOOST_CHECK(fabs(mass - initialMass) / mass < 0.05);
	// NOTE: mass is not conserved exactly for explicit euler (but up to 5 percent)

	cout << "done." << endl;
} /* SEDGMinLee_streaming_test */

/*
 BOOST_AUTO_TEST_CASE(SEDGMinLee_unique_dofs_test) {
 // Test if the unique relation between degrees of freedom and quadrature nodes is OK

 cout << "SEDGMinLee_unique_dofs_test..." << endl;
 // Create problem
 size_t refinementLevel = 1;
 size_t fe_order = 3;
 PeriodicTestDomain2D periodic(refinementLevel);
 SEDGMinLee<2> streaming(periodic.getTriangulation(),
 periodic.getBoundaries(), fe_order,
 make_shared<D2Q9IncompressibleModel>());

 // Test if fi(q) == 1 for the given pairs
 dealii::MappingQ1<2> mapping;
 /// integration on gauss lobatto nodes
 shared_ptr<dealii::QGaussLobatto<2> > quadrature = make_shared<
 dealii::QGaussLobatto<2> >(streaming.getOrderOfFiniteElement());
 /// integration on boundary (with gauß lobatto nodes)
 shared_ptr<dealii::QGaussLobatto<1> > faceQuadrature = make_shared<
 dealii::QGaussLobatto<1> >(streaming.getOrderOfFiniteElement());
 shared_ptr<dealii::FE_DGQArbitraryNodes<2> > fe = make_shared<
 dealii::FE_DGQArbitraryNodes<2> >(
 dealii::QGaussLobatto<1>(streaming.getOrderOfFiniteElement()));

 dealii::FEValues<2> feCellValues(mapping, *fe, *quadrature,
 dealii::update_values);
 dealii::FEFaceValues<2> feFaceValues(mapping, *fe, *faceQuadrature,
 dealii::update_values);

 // for cell q points
 feCellValues.reinit(streaming.getDoFHandler()->begin_active());
 std::map<size_t, size_t> dofToQ = streaming.getCelldofToQIndex();
 for (size_t i = 0; i < streaming.getFe()->n_dofs_per_cell(); i++) {
 BOOST_CHECK(fabs(feCellValues.shape_value(i, dofToQ[i]) - 1) < 1e-10);
 }
 // for face q points
 vector<std::map<size_t, size_t> > dofToFaceQ =
 streaming.getFacedofToQIndex();
 for (size_t i = 0; i < 4; i++) {
 feFaceValues.reinit(streaming.getDoFHandler()->begin_active(), i);
 for (size_t j = 0; j < streaming.getFe()->n_dofs_per_cell(); j++) {
 // if key is in map
 if (dofToFaceQ.at(i).count(j) == 1) {
 BOOST_CHECK(
 fabs(
 feFaceValues.shape_value(j, dofToFaceQ.at(i)[j])
 - 1) < 1e-10);
 }
 }
 }

 cout << "done" << endl;
 }  SEDGMinLee_unique_dofs_test */

BOOST_AUTO_TEST_CASE(SEDGMinLee_RKstreaming_test) {
	cout << "SEDGMinLee_RKstreaming_test..." << endl;
/// THIS WILL SHOW, if streaming is really correct

	/// relaxationParamter and velocity have no impact; just needed for construction
	size_t refinementLevel = 3;
	size_t fe_order = 4;
	bool useCentralFlux = false;
	PeriodicTestDomain2D periodic(refinementLevel);

	SEDGMinLee<2> streaming(periodic.getTriangulation(),
			periodic.getBoundaries(), fe_order,
			make_shared<D2Q9IncompressibleModel>(), "", useCentralFlux);
	const distributed_sparse_block_matrix& matrices =
			streaming.getSystemMatrix();

	// choose time step and number of time steps so dx = dt and the bump passes the domain one time (five times)
	const double timeStep = 1. / (pow(2, refinementLevel) * (fe_order - 1));
#ifdef RK5_OUT
	const size_t numberOfTimeSteps = size_t(5. / timeStep);
#else
	const size_t numberOfTimeSteps = 1. / timeStep;
#endif
	const size_t timeStepsUntilHalfThroughDomain = 0.5 / timeStep;

	// Initialize all particle distribution functions with 1, one corner element with 2
	distributed_vector f(streaming.getNumberOfDoFs());
	double mass = 0.0;
	double initialMass = 0.0;
	for (size_t i = 0; i < f.size(); i++) {
		f(i) = 1.0;
		initialMass += 1;
	}
	// create smooth initial conditions
	vector<dealii::Point<2> > supportPoints(
			streaming.getDoFHandler()->n_dofs());
	dealii::DoFTools::map_dofs_to_support_points(streaming.getMapping(),
			*streaming.getDoFHandler(), supportPoints);
	dealii::Point<2> midPoint(0.25, 0.25);
	const double PI = 3.14;
	for (size_t i = 0; i < streaming.getDoFHandler()->n_dofs(); i++) {
		double distance = supportPoints.at(i).distance(midPoint);
		if (distance <= 0.25) {
			f(i) = 1 + 0.1 * cos(PI / 0.5 * distance);
			initialMass += 0.1 * cos(PI / 0.5 * distance);
		}
	}
	// CHECK MASS CONSERVATION
	// RK5 for streaming in direction (-1,0)
	distributed_sparse_matrix advectionMatrix;
	advectionMatrix.reinit(streaming.getSparsityPattern(2));
	advectionMatrix.copy_from(matrices.block(2,2));
	RungeKutta5LowStorage<distributed_sparse_matrix, distributed_vector> RK5(timeStep, f.size());
#ifdef RK5_OUT
	// Make results dir
	std::stringstream dirname;
	dirname << getenv("NATRIUM_HOME") << "/rk5_periodicstreaming";
	if (mkdir(dirname.str().c_str(), 0777) == -1) {
		if (errno != EEXIST) {
			cerr << "Fehler in mkdir: " << strerror(errno) << endl;
		}
	}
#endif
	for (size_t i = 0; i < numberOfTimeSteps; i++) {
#ifdef RK5_OUT
		// output
		std::stringstream str;
		str <<dirname.str() << "/t_" << i << ".dat";
		std::string filename = str.str();
		std::ofstream gnuplot_output(filename.c_str());
		dealii::DataOut<2> data_out;
		data_out.attach_dof_handler(*streaming.getDoFHandler());
		data_out.add_data_vector(f, "f");
		data_out.build_patches();
		data_out.write_gnuplot(gnuplot_output);
#endif
		//stream
		// 1 time steps at once
		const size_t numberOfTimeStepsAtOnce = 1;
		for (size_t j = 0; j < numberOfTimeStepsAtOnce; j++) {
			RK5.step(f, advectionMatrix);
		}

		if (i == timeStepsUntilHalfThroughDomain) {
			// check if bump has passed half of the domain
			dealii::Point<2> expectedBumpMid(0.75, 0.25);
			double maximum = 0.0;
			size_t index = 0;
			for (size_t i = 0; i < f.size(); i++) {
				if (f(i) > maximum) {
					maximum = f(i);
					index = i;
				}
			}
			BOOST_CHECK(fabs(maximum - 1.0) < 0.15);
			double distanceToExpectedBumpMid = supportPoints.at(index).distance(
					expectedBumpMid);
			BOOST_CHECK(distanceToExpectedBumpMid < 0.15);
		}
	}

#ifdef RK5_OUT
	std::stringstream filename;
	filename << getenv("NATRIUM_HOME") << "/rk5_periodicstreaming/matrix.dat";
	std::ofstream matrixfile(filename.str());
	matrices.block(3,3).print_formatted(matrixfile);
#endif

	// check if mass was conserved
	mass = 0.0;
	for (size_t j = 0; j < f.size(); j++) {
		mass += f(j);
	}
	BOOST_CHECK(fabs(mass - initialMass) / initialMass < 1e-3);

	// Check if bump is again at initial position
	double maximum = 0.0;
	size_t index = 0;
	for (size_t i = 0; i < f.size(); i++) {
		if (f(i) > maximum) {
			maximum = f(i);
			index = i;
		}
	}
	BOOST_CHECK(fabs(maximum - 1.0) < 0.15);
	double distToInitialMidPoint = supportPoints.at(index).distance(midPoint);
	BOOST_CHECK(distToInitialMidPoint < 0.15);

	cout << "done." << endl;
} /* SEDGMinLee_RKstreaming_test */


BOOST_AUTO_TEST_CASE(SEDGMinLee_SaveAndLoadCheckpoints_test){
	cout << "SEDGMinLee_SaveAndLoadCheckpoints_test..." << endl;

	std::stringstream directory;
	directory << getenv("NATRIUM_HOME") << "/test-restart";

	// create streaming object
	size_t refinementLevel = 3;
	size_t fe_order = 2;
	bool useCentralFlux = false;
	PeriodicTestDomain2D periodic(refinementLevel);
	SEDGMinLee<2> streaming(periodic.getTriangulation(),
			periodic.getBoundaries(), fe_order,
			make_shared<D2Q9IncompressibleModel>(), "", useCentralFlux);
	streaming.saveCheckpoint(directory.str());

	/////// SANITY TEST //////////
	BOOST_CHECK_NO_THROW(SEDGMinLee<2>(periodic.getTriangulation(),
			periodic.getBoundaries(), fe_order,
			make_shared<D2Q9IncompressibleModel>(), "", useCentralFlux));

	/////// FAILURE TEST ////////
	PeriodicTestDomain2D periodic2(4);
	BOOST_CHECK_THROW(SEDGMinLee<2>(periodic2.getTriangulation(),
			periodic2.getBoundaries(), fe_order,
			make_shared<D2Q9IncompressibleModel>(), directory.str(), useCentralFlux), AdvectionSolverException);

	BOOST_CHECK_THROW(SEDGMinLee<2>(periodic.getTriangulation(),
			periodic.getBoundaries(), fe_order,
			make_shared<D2Q9IncompressibleModel>(), "gubaguba", useCentralFlux), AdvectionSolverException);

	BOOST_CHECK_THROW(SEDGMinLee<2>(periodic.getTriangulation(),
			periodic.getBoundaries(), fe_order+1,
			make_shared<D2Q9IncompressibleModel>(), directory.str(), useCentralFlux), AdvectionSolverException);

	BOOST_CHECK_THROW(SEDGMinLee<2>(periodic.getTriangulation(),
			periodic.getBoundaries(), fe_order,
			make_shared<D2Q9IncompressibleModel>(), directory.str(), not useCentralFlux), AdvectionSolverException);


	cout << "done" << endl;
} /*SEDGMinLee_SaveAndLoadCheckpoints_test*/

BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
