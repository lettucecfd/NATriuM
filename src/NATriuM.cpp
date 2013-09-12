//============================================================================
// Name        : NATriuM.cpp
// Author      : Andreas Kraemer, andreas.kraemer@h-brs.de
// Version     : 0
// Copyright   : University of Applied Sciences Sankt Augustin, Germany
// Description : Numerics and Algorithms for Tribology
//============================================================================

#include <iostream>
#include <vector>

#include "deal.II/grid/tria.h"

#include "boltzmannmodels/D2Q9IncompressibleModel.h"
#include "solver/CFDSolver.h"

using dealii::Triangulation;
using natrium::D2Q9Incompressible;
using natrium::CFDSolver;
using namespace std;


size_t main() {

	cout << "Starting NATriuM." << endl;

	//Triangulation<2> tri;

	vector<natrium::float_t> v(2,0.0);
	D2Q9IncompressibleModel boltzmannModel;


	//CFDSolver<2>* solver = new CFDSolver<2>();



	//float_t rho = 1;
	//cout << boltzmannModel->getEquilibriumDistribution(0,v,rho) << endl;

	cout << "NATriuM terminated." << endl;

}
