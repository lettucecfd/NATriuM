/**
 * @file Logging.h
 * @short Definition of logging output streams
 * @date 19.02.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "Logging.h"

#include <fstream>
#include <sstream>

#include "boost/make_shared.hpp"

namespace natrium {

shared_ptr<Logging> Logging::m_LOGGER;

Logging& LOG(LogLevel level){
	if (0 == Logging::m_LOGGER){
		Logging::m_LOGGER = make_shared<Logging>();
	}
	return (*Logging::m_LOGGER)(level);
}
Logging& LOGGER(){
	if (0 == Logging::m_LOGGER){
		Logging::m_LOGGER = make_shared<Logging>();
	}
	return *Logging::m_LOGGER;
}

const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);

    return buf;
}

} /* namespace natrium */
