#ifndef SRC_METRIC_H_
#define SRC_METRIC_H_

#include <iostream>
#include <stdexcept>
#include <limits>
#include <vector>
#include <string>
#include "aux.h"
#include "GzipWriter.h"
#include "Histogram.h"

using namespace std;

class Metric {
private:
	double* data;
	unsigned int allocated;
	unsigned int n;

	string name;
	string description;

public:
	Metric(const string& name, const string& description);
	Metric(Metric &&metric);

	void add(double value) throw (runtime_error);
	unsigned int get_length();
	double get(unsigned int i);
	void write_histogram(GzipWriter& writer, unsigned int n_bins) throw (runtime_error);
	static void write_histogram(GzipWriter& writer, vector<Metric> &metrics, unsigned int n_bins) throw (runtime_error);

	virtual ~Metric();
};

#endif
