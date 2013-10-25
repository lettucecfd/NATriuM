/**
 * @file step-1.cpp
 * @short First tutorial:  Couette Flow in 2D
 * @date 24.10.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "CouetteFlow2D.h"

#include "utilities/BasicNames.h"


using namespace natrium;

// Main function
int main() {

	cout << "Starting NATriuM step-1..." << endl;

	double relaxationParameter = 0.7;
	double topPlateVelocity = 0.05;

	CouetteFlow2D couetteFlow(relaxationParameter, topPlateVelocity);

	cout << "NATriuM step-1 terminated." << endl;

	return 0;
}
