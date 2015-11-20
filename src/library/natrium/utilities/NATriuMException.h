/*
 * NATriuMException.h
 *
 *  Created on: Feb 26, 2015
 *      Author: kraemer
 */

#ifndef NATRIUMEXCEPTION_H_
#define NATRIUMEXCEPTION_H_

#include <exception>
#include <mpi.h>

#include "BasicNames.h"
#include "../utilities/Logging.h"

namespace natrium {

/**
 * @short Exception class for CFDSolver
 */
class NATriuMException: public std::exception {
private:
	std::string message;
public:
	NATriuMException(const char *msg) :
			message(msg) {
		LOG(ERROR) << "ERROR in NATriuM. Error message:" << msg << endl;
	}
	NATriuMException(const std::string& msg) :
			message(msg) {
	}
	~NATriuMException() throw () {
	}
	virtual const char *what() const throw () {
		return this->message.c_str();
	}
};

inline void natrium_errorexit(const char* msg) {

	LOG(ERROR) << "---------------------------------------------" << endl;
	LOG(ERROR) << "An error occurred in your NATriuM simulation." << endl;
	if (dealii::Utilities::MPI::job_supports_mpi())
		LOG(ERROR) << "Killing all MPI processes." << endl;
	LOG(ERROR) << "The error message:" << endl << msg << endl;
	LOG(ERROR) << "---------------------------------------------" << endl;
	LOG(ERROR) << "Simulation failed" << endl;

	pout << "Error:" << msg << endl;

	// Kill all MPI jobs
	if (dealii::Utilities::MPI::job_supports_mpi())
		MPI_Abort(MPI_COMM_WORLD, 1);
}

} /* namespace  natrium */

#endif /* NATRIUMEXCEPTION_H_ */
