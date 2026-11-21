// Frequently Used Functions

#pragma once

#include "base/Global.hpp"
#include "base/Array.hpp"

template <typename type>
void print(type value) {
	std::cout << value << std::endl;
}

template <typename type>
void print(std::vector<type>& data) {
	for (size_t i = 0; i < data.size(); i++) {
		std::cout << data[i] << " ";
	}
	std::cout << std::endl;
}

template <typename type>
type find_max(const type& a, const type& b) {
	if (std::isnan(a)) { return b; }
	else if (std::isnan(b)) { return a; }
	else if (a > b) { return a; }
	else { return b; }
}

template <typename type>
type find_max(Array<type>& data, label dim = -1) {
	type max = 0;
	if (dim == -1) {
		max = data(0);
		for (label i = 0; i < data.size; i++) {
			if (data(i) > max) { max = data(i); }
		}
	}
	else {
		max = data(0, dim);
		for (label i = 0; i < data.NX; i++) {
			if (data(i, dim) > max) { max = data(i, dim); }
		}
	}
	return max;
}

template <typename type>
type find_min(const type& a, const type& b) {
	if (std::isnan(a)) { return b; }
	else if (std::isnan(b)) { return a; }
	else if (a < b) { return a; }
	else { return b; }
}

bool getline_clean(std::istream& file, string& line);

scalar fabs_max(scalar& a, scalar& b);
