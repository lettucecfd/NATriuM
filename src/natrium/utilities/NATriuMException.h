/*
 * NATriuMException.h
 *
 *  Created on: Feb 26, 2015
 *      Author: kraemer
 */

#ifndef NATRIUMEXCEPTION_H_
#define NATRIUMEXCEPTION_H_

#include <exception>

#include "BasicNames.h"
#include "../utilities/Logging.h"

namespace  natrium {

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
	const char *what() const throw () {
		return this->message.c_str();
	}
};

} /* namespace  natrium */

#endif /* NATRIUMEXCEPTION_H_ */
