#include <string>
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

// code has to be called from a git directory to make this function work
std::string getGitSha() {
		return exec("git show --abbrev-commit 2> NULL | head -1 | awk '{print $2;}'");
}

std::string getUserName() {
	return exec("whoami");
}

std::string getHostName() {
	return exec("hostname");
}

} /* namespace Info */

} /* namespace natrium*/
