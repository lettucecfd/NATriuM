/**
 * @file SolverConfiguration.h
 * @short Class that stores the configuration for a CFD simulation based on the Discrete Boltzmann Equation (DBE).
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef SOLVERCONFIGURATION_H_
#define SOLVERCONFIGURATION_H_

namespace natrium {

/** @short Class that stores the configuration for a CFD simulation based on the Discrete Boltzmann Equation (DBE).
 */
class SolverConfiguration {
private:
public:

	/// constructor
	SolverConfiguration();

	/// destructor
	virtual ~SolverConfiguration();
};

} /* namespace natrium */
#endif /* SOLVERCONFIGURATION_H_ */
