/**
 * @file DataMinLee2011_test.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "streamingdata/DataMinLee2011.h"

#include <fstream>
#include "complex.h"

#include "boost/test/unit_test.hpp"

#include "deal.II/numerics/data_out.h"
#include "deal.II/dofs/dof_tools.h"
#include "deal.II/lac/lapack_full_matrix.h"

#include "boltzmannmodels/D2Q9IncompressibleModel.h"

#include "timeintegration/RungeKutta5LowStorage.h"

#include "PeriodicTestDomain2D.h"

#include "utilities/BasicNames.h"

///////////////////////////////
///////////////////////////////
// FLAGS TO SWITCH OUTPUT ON //
//#define CREATE_DATA_FILES // This takes a lot of time (~1h)
//#define PRINT_SYSTEM_MATRIX
//#define EULER_OUT
#define RK5_OUT
///////////////////////////////
///////////////////////////////

namespace natrium {

BOOST_AUTO_TEST_SUITE(DataMinLee2011_test)

BOOST_AUTO_TEST_CASE(DataMinLee2011_Construction_test) {
	cout << "DataMinLee2011_Construction_test..." << endl;

	size_t fe_order = 2;
	size_t refinementLevel = 2;
	PeriodicTestDomain2D periodic(refinementLevel);
	BOOST_CHECK_NO_THROW(
			DataMinLee2011<2> streaming(periodic.getTriangulation(), periodic.getBoundaries(), fe_order, make_shared<D2Q9IncompressibleModel>()));

	cout << "done." << endl;
} /* DataMinLee2011_Construction_test */

BOOST_AUTO_TEST_CASE(DataMinLee2011_systemMatrix_test) {
	cout << "DataMinLee2011_systemMatrix_test..." << endl;

	bool useCentralFlux = true;
	do { // useCentralFlux = true/false
#ifdef CREATE_DATA_FILES
	for (size_t fe_order = 2; fe_order < 10; fe_order++) {
		// create files for max eigenvalue diagram
		std::stringstream maxEigenvalueFilename;
		maxEigenvalueFilename << "../results/maxeigenvalue_advection/";
		if (useCentralFlux) {
			maxEigenvalueFilename << "maxeigenvalue_cflux_";
		} else {
			maxEigenvalueFilename << "maxeigenvalue_laxflux_";
		}
		maxEigenvalueFilename << "feorder_" << fe_order << ".dat";
		cout << maxEigenvalueFilename.str().c_str() << "..." << endl;
		std::ofstream normOut(maxEigenvalueFilename.str().c_str());
#else
		for (size_t fe_order = 2; fe_order <3; fe_order ++) {
#endif
			for (size_t refinementLevel = 1; refinementLevel < 4;
					refinementLevel++) {
				PeriodicTestDomain2D periodic(refinementLevel);
				DataMinLee2011<2> streaming(periodic.getTriangulation(),
						periodic.getBoundaries(), fe_order,
						make_shared<D2Q9IncompressibleModel>(), useCentralFlux);
				const vector<distributed_sparse_matrix>& matrices =
						streaming.getSystemMatrix();
				BOOST_CHECK(matrices.size() == D2Q9IncompressibleModel::Q);
#ifdef PRINT_SYSTEM_MATRIX
				for (size_t i = 0; i < 2; i++) { //D2Q9IncompressibleModel::Q; i++){
					cout << "Matrix " << i << ": " << endl;
					matrices.at(i).print_formatted(cout);
				}
#endif
				//BOOST_CHECK()
				//singular value decomposition
				dealii::LAPACKFullMatrix<double> systemMatrix(
						matrices.at(0).n());
				for (size_t i = 0; i < D2Q9IncompressibleModel::Q; i++) {
#ifdef CREATE_DATA_FILES
					// create files for spectrum plots
					std::stringstream filename;
					filename << "../results/eigenvalues_advection/";
					if (useCentralFlux) {
						filename << "cflux_";
					} else {
						filename << "laxflux_";
					}
					filename << "refine_" << refinementLevel << "_order_"
					<< fe_order << "_matrix_" << i << "_spectrum.dat";
					cout << filename.str().c_str() << "..." << endl;
					std::ofstream spectrumOut(filename.str().c_str());
#endif
					// compute eigenvalues
					systemMatrix = 0;
					systemMatrix.copy_from(matrices.at(i));
					systemMatrix.compute_eigenvalues();
					// check eigenvalues (the real part must vanish for central flux, the absolute value must be in the same order of magnitude as delta x)
					double max_eigenvalue = 0;
					for (size_t j = 0; j < systemMatrix.n_cols(); j++) {
						std::
						complex<double> eigenValue;
						eigenValue = systemMatrix.eigenvalue(j);
						BOOST_CHECK(eigenValue.real() < 0.1);
						BOOST_CHECK(eigenValue.imag() < 1.0);
						if (useCentralFlux) {
							BOOST_CHECK(fabs(eigenValue.real()) < 1e-15);
						}
						if (abs(eigenValue) > max_eigenvalue) {
							max_eigenvalue = abs(eigenValue);
						}
					}
					double delta_x = 1. / pow(2, refinementLevel + 1);
					if (i != 0) {
						BOOST_CHECK(max_eigenvalue > 0.5 * delta_x);
						BOOST_CHECK(max_eigenvalue < 2.5 * delta_x);
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
		useCentralFlux = !useCentralFlux;
	} while (!useCentralFlux);

	cout << "done." << endl;
}
/* DataMinLee2011_systemMatrix_test */

BOOST_AUTO_TEST_CASE(DataMinLee2011_steadyStreaming_test) {
	cout << "DataMinLee2011_steadyStreaming_test..." << endl;

	size_t refinementLevel = 4;
	size_t fe_order = 4;
	PeriodicTestDomain2D periodic(refinementLevel);
	DataMinLee2011<2> streaming(periodic.getTriangulation(),
			periodic.getBoundaries(), fe_order,
			make_shared<D2Q9IncompressibleModel>());
	const vector<distributed_sparse_matrix>& matrices =
			streaming.getSystemMatrix();
	const double timeStep = 0.05;
	const size_t numberOfTimeSteps = 10;

	// Initialize all particle distribution functions with 1
	distributed_vector f(streaming.getSystemMatrix().at(0).n());
	for (size_t i = 0; i < f.size(); i++) {
		f(i) = 1.0;
	}

	// CHECK MASS CONSERVATION AND CONSERVATION OF STEADY STATE
	double initialMass = f.size();
	// Explicit euler for streaming in direction (1,0)
	distributed_sparse_matrix advectionMatrix;
	advectionMatrix.reinit(streaming.getSparsityPattern());
	advectionMatrix.copy_from(matrices.at(3));
	advectionMatrix *= timeStep;
	distributed_vector f_tmp(streaming.getSystemMatrix().at(0).n());
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
} /* DataMinLee2011_steadyStreaming_test */

BOOST_AUTO_TEST_CASE(DataMinLee2011_streaming_test) {
	cout << "DataMinLee2011_streaming_test..." << endl;

	size_t refinementLevel = 4;
	size_t fe_order = 4;
	PeriodicTestDomain2D periodic(refinementLevel);
	DataMinLee2011<2> streaming(periodic.getTriangulation(),
			periodic.getBoundaries(), fe_order,
			make_shared<D2Q9IncompressibleModel>());
	const vector<distributed_sparse_matrix>& matrices =
			streaming.getSystemMatrix();
	const double timeStep = 0.5;
	const size_t numberOfTimeSteps = 50;

	// Initialize all particle distribution functions with 1, one corner element with 2
	distributed_vector f(streaming.getSystemMatrix().at(0).n());
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
	advectionMatrix.reinit(streaming.getSparsityPattern());
	advectionMatrix.copy_from(matrices.at(3));
	advectionMatrix *= timeStep;
	distributed_vector f_tmp(streaming.getSystemMatrix().at(0).n());
	for (size_t i = 0; i < numberOfTimeSteps; i++) {
#ifdef EULER_OUT
		// output
		std::stringstream str;
		str << "../results/expliciteuler_periodicstreaming/t_" << i << ".dat";
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
<<<<<<< HEAD
		// check if mass is conserved
=======

		double mass = 0.0;
		for (size_t j = 0; j < f.size(); j++) {
			mass += f(j);
		}

		BOOST_CHECK(fabs(mass - initialMass) / mass < 1e-2);
		// NOTE: mass is not conserved exactly for explicit euler (but up to 5 percent)

>>>>>>> 889e86e40cdea2e1a4b0fa8848567c004e24cdbf
	}
	double mass = 0.0;
	for (size_t j = 0; j < f.size(); j++) {
		mass += f(j);
	}
	BOOST_CHECK(fabs(mass - initialMass) / mass < 0.05);
	// NOTE: mass is not conserved exactly for explicit euler (but up to 5 percent)

	cout << "done." << endl;
} /* DataMinLee2011_streaming_test */

BOOST_AUTO_TEST_CASE(DataMinLee2011_RKstreaming_test) {
	cout << "DataMinLee2011_RKstreaming_test..." << endl;
	/// THIS WILL SHOW, if streaming is really correct

	/// relaxationParamter and velocity have no impact; just needed for construction
	size_t refinementLevel = 4;
	size_t fe_order = 2;
	bool useCentralFlux = false;
	PeriodicTestDomain2D periodic(refinementLevel);
	DataMinLee2011<2> streaming(periodic.getTriangulation(),
			periodic.getBoundaries(), fe_order,
			make_shared<D2Q9IncompressibleModel>(), useCentralFlux);
	const vector<distributed_sparse_matrix>& matrices =
			streaming.getSystemMatrix();
	const double timeStep = 1;
<<<<<<< HEAD
#ifdef RK5_OUT
	// The videoplot.sh file take by default 500 files
=======
>>>>>>> 889e86e40cdea2e1a4b0fa8848567c004e24cdbf
	const size_t numberOfTimeSteps = 500;
#else
	const size_t numberOfTimeSteps = 50;
#endif
	// Initialize all particle distribution functions with 1, one corner element with 2
	distributed_vector f(streaming.getSystemMatrix().at(0).n());
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
	// RK5 for streaming in direction (1,0)
	distributed_sparse_matrix advectionMatrix;
	advectionMatrix.reinit(streaming.getSparsityPattern());
	advectionMatrix.copy_from(matrices.at(3));
	RungeKutta5LowStorage RK5(timeStep, f.size());
	for (size_t i = 0; i < numberOfTimeSteps; i++) {
#ifdef RK5_OUT
		// output
		std::stringstream str;
		str << "../results/rk5_periodicstreaming/t_" << i << ".dat";
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
		// check if mass is conserved
		double mass = 0.0;
		for (size_t j = 0; j < f.size(); j++) {
			mass += f(j);
		}
		BOOST_CHECK(fabs(mass - initialMass) / mass < 1e-5);
	}

	cout << "done." << endl;
} /* DataMinLee2011_RKstreaming_test */

BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
