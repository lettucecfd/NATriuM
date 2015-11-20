
#ifndef NATRIUM_INFO_H_
#define NATRIUM_INFO_H_

#include <boost/filesystem.hpp>

#include "BasicNames.h"

namespace natrium {

namespace Info {

/**
 * @short execute a command on the command line
 * @param[in] cmd A string that contains the command.
 * @return A string that contains the response of the command line
 */
std::string exec(char const * cmd) ;

/**
 * @short get NATRIUM_DIR environment variable
 * @note The environement variables NATRIUM_DIR and NATRIUM_HOME have to be set to make NATriuM work.
 */
boost::filesystem::path get_natrium_dir();

/**
 * @short Return identifier for the current git commit.
 */
std::string getGitSha() ;

std::string getUserName();

std::string getHostName() ;

} /* namespace Info */

} /* namespace natrium*/

#endif /* ifndef NATRIUM_INFO_H */
