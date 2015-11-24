/*
 * Timing.h
 *
 *  Created on: 24.11.2015
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_UTILITIES_TIMING_H_
#define LIBRARY_NATRIUM_UTILITIES_TIMING_H_

#include <mpi.h>

#include "deal.II/base/timer.h"

#include "BasicNames.h"

using dealii::TimerOutput;

namespace natrium {

/**
 * @short Kind of a singleton for Runtime measurements
 */
class Timing {
private:
	/// outstream to print not to cout but to Logfile
	static shared_ptr<std::ostringstream> m_outStream;
	/// TimerOutput is deal.II's class for time measurements
	static shared_ptr<TimerOutput> m_timer;
	/// private constructor
	Timing() {
	}
public:
	virtual ~Timing() {

	}
	/**
	 * @short public constructor; initializes static instance
	 */
	static void makeInstance();
	/**
	 * @short static function to get the global timer
	 */
	static TimerOutput& getTimer();
	/**
	 * @short Static function to get the global outstream.
	 * 		  The outstream is printed to via print_summary.
	 * 		  It can be used to redirect the output, e.g. to the Logfile.
	 *
	 */
	static std::ostringstream& getOutStream();
};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_UTILITIES_TIMING_H_ */
