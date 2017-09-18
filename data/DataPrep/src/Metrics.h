#ifndef SRC_METRICS_H_
#define SRC_METRICS_H_

#include <iostream>
#include <stdexcept>
#include <limits>
#include <vector>
#include <string>
#include "aux.h"
#include "GzipWriter.h"
#include "Histogram.h"

using namespace std;

class Metrics {
private:
	double* data;
	unsigned int allocated;
	unsigned int n;

	string name;

public:
	Metrics(const string& name);

	void add(double value) throw (runtime_error);
	unsigned int get_length();
	double get(unsigned int i);
	void write_histogram(GzipWriter& writer, unsigned int n_bins) throw (runtime_error);

	virtual ~Metrics();
};

#endif
