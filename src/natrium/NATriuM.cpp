//============================================================================
// Name        : NATriuM.cpp
// Author      : Andreas Kraemer, andreas.kraemer@h-brs.de
// Version     : 0
// Copyright   : University of Applied Sciences Sankt Augustin, Germany
// Description : Numerics and Algorithms for Tribology
//============================================================================

#include "deal.II/grid/tria.h"

#include "boltzmannmodels/D2Q9IncompressibleModel.h"
#include "solver/CFDSolver.h"

#include "utilities/BasicNames.h"

using namespace natrium;
//using dealii::Triangulation;


int main() {

	cout << "Starting NATriuM." << endl;

	//Triangulation<2> tri;

	vector<double> v(2,0.0);
	D2Q9IncompressibleModel boltzmannModel;


	//CFDSolver<2>* solver = new CFDSolver<2>();



	//double rho = 1;
	//cout << boltzmannModel->getEquilibriumDistribution(0,v,rho) << endl;

	cout << "NATriuM terminated." << endl;

	return 0;

}
