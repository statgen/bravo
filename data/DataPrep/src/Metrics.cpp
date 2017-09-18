#include "Metrics.h"

Metrics::Metrics(const string& name): data(nullptr), allocated(0u), n(0u), name(name) {

}

Metrics::~Metrics() {
	if (data != nullptr) {
		free(data);
		data = nullptr;
	}
}

void Metrics::add(double value) throw (runtime_error) {
	if (n >= allocated) {
		double* reallocated_data = nullptr;
		allocated += 10000000u;
		if ((reallocated_data = (double*)realloc(data, allocated * sizeof(double))) == nullptr) {
			throw runtime_error("Error while allocating memory!");
		}
		data = reallocated_data;
		reallocated_data = nullptr;
	}
	data[n++] = value;
}

unsigned int Metrics::get_length() {
	return n;
}

double Metrics::get(unsigned int i) {
	return data[i];
}

void Metrics::write_histogram(GzipWriter& writer, unsigned int n_bins) throw (runtime_error) {
	double min = numeric_limits<double>::max();
	double max = numeric_limits<double>::min();

	for (unsigned int i = 0u; i < n; ++i) {
		if (aux::fcmp(min, data[i], 0.00000001) > 0) {
			min = data[i];
		}
		if (aux::fcmp(max, data[i], 0.00000001) < 0) {
			max = data[i];
		}
	}

	vector<double> hist_borders;

	if ((n_bins > 1) && (aux::fcmp(min, max, 0.00000001) != 0)) {
		double step = (max - min) / ((double) n_bins);
		double border = min;
		while (aux::fcmp(border, max, 0.00000001) <= 0) {
			hist_borders.push_back(border);
			border += step;
		}
	} else {
		hist_borders.push_back(min);
		hist_borders.push_back(max);
	}

	Histogram histogram(hist_borders);

	for (unsigned int i = 0u; i < n; ++i) {
		histogram.add(data[i]);
	}

	writer.write("{");
	writer.write("\"metric\":%s,", name.c_str());
	writer.write("\"n\":%d,", n);
	writer.write("\"min\":%g,\"max\":%g,", min, max);
	writer.write("\"hist\":[");
	vector<pair<double, unsigned int>> bins = histogram.get_bins();

	vector<pair<double, unsigned int>>::iterator bins_it = bins.begin();
	writer.write("%d", bins_it->second);
	while(++bins_it != bins.end() - 1) {
		writer.write(",%d", bins_it->second);
	}
	writer.write("],");
	writer.write("\"mids\":[");
	bins_it = bins.begin() + 1;
	writer.write("%g", (bins_it->first - (bins_it->first - (bins_it - 1)->first) / 2.0));
	while (++bins_it != bins.end()) {
		writer.write(",%g", (bins_it->first - (bins_it->first - (bins_it - 1)->first) / 2.0));
	}
	writer.write("]");
	writer.write("}");
}

