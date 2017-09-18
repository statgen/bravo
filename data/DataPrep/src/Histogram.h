#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_

#include <iostream>
#include <vector>
#include <set>
#include <utility>
#include <algorithm>
#include <sstream>
#include <stdexcept>
#include "aux.h"

using namespace std;

class Histogram {
private:
	vector<pair<double, unsigned int>> bins;
	vector<double> borders;

	pair<double, unsigned int> last_added_value;

	char* text;

public:
	Histogram(const vector<double>& borders);
	Histogram(Histogram &&histogram);
	virtual ~Histogram() noexcept;

	void add(double value) noexcept;
	vector<pair<double, unsigned int>> get_bins() noexcept;
	void clear() noexcept;
	const char* get_text() throw (runtime_error);

};

#endif
