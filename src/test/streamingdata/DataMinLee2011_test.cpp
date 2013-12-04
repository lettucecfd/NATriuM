/**
 * @file DataMinLee2011_test.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "streamingdata/DataMinLee2011.h"

#include <fstream>
//#include <strstream>

#include "boost/test/unit_test.hpp"

#include "deal.II/numerics/data_out.h"

#include "boltzmannmodels/D2Q9IncompressibleModel.h"

#include "PeriodicFlow2D.h"

#include "utilities/BasicNames.h"


namespace natrium {

BOOST_AUTO_TEST_SUITE(DataMinLee2011_test)

BOOST_AUTO_TEST_CASE(DataMinLee2011_Construction_test) {
	cout << "DataMinLee2011_Construction_test..." << endl;

	double relaxationParameter = 0.7;
	numeric_vector velocity(2);
	velocity(0) = 0.05;
	velocity(1) = 0.01;
	size_t fe_order = 2;
	PeriodicFlow2D periodic(relaxationParameter, velocity);
	BOOST_CHECK_NO_THROW(DataMinLee2011<2> streaming(periodic.getTriangulation(), periodic.getBoundaries(), fe_order, make_shared<D2Q9IncompressibleModel>()));

	cout << "done." << endl;
} /* DataMinLee2011_Construction_test */

BOOST_AUTO_TEST_CASE(DataMinLee2011_systemMatrix_test) {
	cout << "DataMinLee2011_systemMatrix_test..." << endl;

	double relaxationParameter = 0.7;
	numeric_vector velocity(2);
	velocity(0) = 0.05;
	velocity(1) = 0.01;
	size_t fe_order = 2;
	PeriodicFlow2D periodic(relaxationParameter, velocity);
	DataMinLee2011<2> streaming(periodic.getTriangulation(), periodic.getBoundaries(), fe_order, make_shared<D2Q9IncompressibleModel>());
	const vector<distributed_sparse_matrix>& matrices = streaming.getSystemMatrix();
	BOOST_CHECK(matrices.size() == D2Q9IncompressibleModel::Q);
	/*for (size_t i = 0; i < 2; i++){//D2Q9IncompressibleModel::Q; i++){
		cout << "Matrix " << i << ": " << endl;
		matrices.at(i).print_formatted(cout);
	}*/

	cout << "done." << endl;
} /* DataMinLee2011_systemMatrix_test */

BOOST_AUTO_TEST_CASE(DataMinLee2011_steadyStreaming_test) {
	cout << "DataMinLee2011_steadyStreaming_test..." << endl;

	double relaxationParameter = 0.7;
	numeric_vector velocity(2);
	velocity(0) = 0.05;
	velocity(1) = 0.01;
	size_t fe_order = 4;
	PeriodicFlow2D periodic(relaxationParameter, velocity);
	DataMinLee2011<2> streaming(periodic.getTriangulation(), periodic.getBoundaries(), fe_order, make_shared<D2Q9IncompressibleModel>());
	const vector<distributed_sparse_matrix>& matrices = streaming.getSystemMatrix();
	const double timeStep = 10;
	const size_t numberOfTimeSteps = 100;

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
	for (size_t i = 0; i < numberOfTimeSteps; i++){
		//stream
		f_tmp = f;
		advectionMatrix.vmult_add(f, f_tmp);
		double mass = 0.0;
		for (size_t j = 0; j < f.size(); j++){
			BOOST_CHECK(abs(f(j)-1.0) < 1e-5);
			mass += f(j);
		}
		BOOST_CHECK(abs(mass-initialMass) < 1e-10);

	}

	cout << "done." << endl;
} /* DataMinLee2011_steadyStreaming_test */

BOOST_AUTO_TEST_CASE(DataMinLee2011_streaming_test) {
	cout << "DataMinLee2011_streaming_test..." << endl;

	double relaxationParameter = 0.7;
	numeric_vector velocity(2);
	velocity(0) = 0.05;
	velocity(1) = 0.01;
	size_t fe_order = 4;
	PeriodicFlow2D periodic(relaxationParameter, velocity);
	DataMinLee2011<2> streaming(periodic.getTriangulation(), periodic.getBoundaries(), fe_order, make_shared<D2Q9IncompressibleModel>());
	const vector<distributed_sparse_matrix>& matrices = streaming.getSystemMatrix();
	const double timeStep = 1.;
	const size_t numberOfTimeSteps = 100;

	// Initialize all particle distribution functions with 1, one corner element with 2
	distributed_vector f(streaming.getSystemMatrix().at(0).n());
	double initialMass = 0.0;
	for (size_t i = 0; i < f.size(); i++) {
			f(i) = 1.0;
			initialMass += 1;
	}
	for (size_t i = 0; i < 9; i++){
		f(i) = 2.0;
		initialMass += 1;
	}
	// CHECK MASS CONSERVATION
	// Explicit euler for streaming in direction (1,0)
	distributed_sparse_matrix advectionMatrix;
	advectionMatrix.reinit(streaming.getSparsityPattern());
	advectionMatrix.copy_from(matrices.at(3));
	advectionMatrix *= timeStep;
	distributed_vector f_tmp(streaming.getSystemMatrix().at(0).n());
	for (size_t i = 0; i < numberOfTimeSteps; i++){

		// output
		std::stringstream str;
		str << "t_" << i << ".gnp";
		std::string filename = str.str();
		std::ofstream gnuplot_output (filename.c_str());
		dealii::DataOut<2> data_out;
		data_out.attach_dof_handler (*streaming.getDoFHandler());
		data_out.add_data_vector (f, "f");
		data_out.build_patches ();
		data_out.write_gnuplot(gnuplot_output);

		//stream
		f_tmp = f;
		advectionMatrix.vmult_add(f, f_tmp);

		double mass = 0.0;
		for (size_t j = 0; j < f.size(); j++){
			mass += f(j);
		}

		BOOST_CHECK(abs(mass-initialMass)/mass < 0.05);
		// NOTE: mass is not conserved exactly for explicit euler (but up to 5 percent)

	}

	cout << "done." << endl;
} /* DataMinLee2011_streaming_test */

BOOST_AUTO_TEST_CASE(DataMinLee2011_RKstreaming_test) {
	cout << "DataMinLee2011_RKstreaming_test..." << endl;
	/// THIS WILL SHOW, if streaming is really correct
	cout << "done." << endl;
} /* DataMinLee2011_RKstreaming_test */


BOOST_AUTO_TEST_SUITE_END()

} /* namespace natrium */
