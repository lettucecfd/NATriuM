
#include <algorithm>
#include <iostream>
#include <stdio.h>

namespace natrium {

namespace Info {

std::string exec(char const * cmd) {
	FILE* pipe = popen(cmd, "r");
	if (!pipe)
		return "ERROR";
	char buffer[128];
	std::string result = "";
	while (!feof(pipe)) {
		if (fgets(buffer, 128, pipe) != NULL)
			result += buffer;
	}
	pclose(pipe);
	result.erase(std::remove(result.begin(), result.end(), '\n'), result.end());
	return result;
}

boost::filesystem::path get_natrium_dir(){
	std::string natrium_dir = getenv("NATRIUM_DIR");
	if (0 == natrium_dir.size()){
		perr << "You have to specify the environment variable NATRIUM_DIR. Error." << endl;
		assert(false);
	}
	boost::filesystem::path result(natrium_dir);
	return result;
}

// code has to be called from a git directory to make this function work
std::string getGitSha() {
		boost::filesystem::path working_dir( boost::filesystem::current_path() );
		boost::filesystem::current_path( get_natrium_dir() );
		std::string result = exec("git show --abbrev-commit 2> /dev/null | head -1 | awk '{print $2;}'");
		boost::filesystem::current_path ( working_dir );
		return result;
}

// code has to be called from a git directory to make this function work
std::string getGitBranch() {
    boost::filesystem::path working_dir( boost::filesystem::current_path() );
    boost::filesystem::current_path( get_natrium_dir() );
    std::string result = exec("git branch | grep '*' | sed 's/^..//'");
    boost::filesystem::current_path ( working_dir );
    return result;
}  //  GIT_BRANCH=`git branch | grep "^\*" | sed 's/^..//'`

std::string getUserName() {
	return exec("whoami");
}

std::string getHostName() {
	return exec("hostname");
}

} /* namespace Info */

} /* namespace natrium*/
