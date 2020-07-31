#ifndef AUX_H_
#define AUX_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <cstdarg>
#include <regex>
#include <algorithm>
#include <map>
#include "bgzf.h"

using namespace std;

namespace aux {

    /* to read samples from list saved in the text file - one sample per line */
    string read_samples(const char* samples_file) throw (runtime_error);

    /* to read samples from tab-delimited (without header) table saved in the text file */
    string read_samples(const char* samples_file, unsigned int col_idx) throw (runtime_error);

    /* to read sample's current name and new name from tab-delimited (without header) table saved in the text file */
    map<string, string> read_sample_names_map(const char* samples_file, unsigned int name_col_idx, unsigned int new_name_col_idx) throw (runtime_error);

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
