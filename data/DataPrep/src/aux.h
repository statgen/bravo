#ifndef AUX_H_
#define AUX_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <cstdarg>
#include "bgzf.h"

using namespace std;

namespace aux {

	string read_samples(const char* samples_file) throw (runtime_error);

	void write(BGZF* f, const char* format, ...) throw (runtime_error);

	inline int fcmp(double x, double y, double epsilon) {
		int max_exponent = 0;
		double delta = 0.0;
		double diff = 0.0;

		frexp(fabs(x) > fabs(y) ? x : y, &max_exponent);
		delta = ldexp(epsilon, max_exponent);

		diff = x - y;

		if (diff > delta) {
			return 1;
		} else if (diff < -delta) {
			return -1;
		} else {
			return 0;
		}
	}

}

#endif
