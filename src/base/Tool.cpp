#include "base/Tool.hpp"

bool getline_clean(std::istream& file, string& line) {  // getline and remove < > and <'> while getline 
	bool success = static_cast<bool>(getline(file, line));
	if (success) {
		line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
		line.erase(std::remove_if(line.begin(), line.end(), [](char c) { return c == '\''; }), line.end());
	}
	return success;
}

scalar fabs_max(scalar& a, scalar& b) {
	if (fabs(a) >= fabs(b)) { return fabs(a); }
	else { return fabs(b); }
}
