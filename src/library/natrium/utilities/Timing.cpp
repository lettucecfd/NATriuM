/*
 * Timing.cpp
 *
 *  Created on: 24.11.2015
 *      Author: akraem3m
 */

#include "Timing.h"

namespace natrium {

boost::shared_ptr<TimerOutput> Timing::m_timer = NULL;
boost::shared_ptr<std::ostringstream> Timing::m_outStream = NULL;

void Timing::makeInstance() {
	if (m_timer == NULL) {
		assert(m_outStream == NULL);
		m_outStream = boost::make_shared<std::ostringstream>();
		std::ostringstream& out = *(m_outStream);
		m_timer = boost::make_shared<TimerOutput>(MPI_COMM_WORLD, out,
				TimerOutput::summary, TimerOutput::wall_times);

	}
}

TimerOutput& Timing::getTimer() {
	makeInstance();
	return *m_timer;
}
/**
 * @short Static function to get the global outstream.
 * 		  The outstream is printed to via print_summary.
 * 		  It can be used to redirect the output, e.g. to the Logfile.
 *
 */
std::ostringstream& Timing::getOutStream() {
	makeInstance();
	return *m_outStream;
}

} /* namespace natrium */
