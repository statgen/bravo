#include "Metric.h"

Metric::Metric(const string& name, const string& description): data(nullptr), allocated(0u), n(0u), name(name), description(description) {

}

Metric::Metric(Metric &&metric) {

	this->data = metric.data;
	metric.data = nullptr;
	this->allocated = metric.allocated;
	this->n = metric.n;
	this->name = std::move(metric.name);
	this->description = std::move(metric.description);
}

Metric::~Metric() {

	if (data != nullptr) {
		free(data);
		data = nullptr;
	}
}

void Metric::add(double value) throw (runtime_error) {
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

unsigned int Metric::get_length() {
	return n;
}

double Metric::get(unsigned int i) {
	return data[i];
}

void Metric::write_histogram(GzipWriter& writer, unsigned int n_bins) throw (runtime_error) {
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
	writer.write("\"description\":\"%s\",", description.c_str());
	writer.write("\"n\":%d,", n);
	writer.write("\"min\":%g,\"max\":%g,", min, max);
	writer.write("\"hist\":[");
	vector<pair<double, unsigned int>> bins = histogram.get_bins();
	vector<pair<double, unsigned int>>::iterator bins_it = bins.begin();
	writer.write("%d", bins_it->second);
	if (bins.size() > 1) {
		while(++bins_it != bins.end() - 1) {
			writer.write(",%d", bins_it->second);
		}
	}
	writer.write("],");
	writer.write("\"mids\":[");
	if (bins.size() > 1) {
		bins_it = bins.begin() + 1;
		writer.write("%g", (bins_it->first - (bins_it->first - (bins_it - 1)->first) / 2.0));
		while (++bins_it != bins.end()) {
			writer.write(",%g", (bins_it->first - (bins_it->first - (bins_it - 1)->first) / 2.0));
		}
	} else {
		writer.write("%g", bins.front());
	}
	writer.write("]");
	writer.write("}");
}

void Metric::write_histogram(GzipWriter& writer, vector<Metric> &metrics, unsigned int n_bins) throw (runtime_error) {
	double min = numeric_limits<double>::max();
	double max = numeric_limits<double>::min();
	unsigned int n = 0u;

	for (auto&& metric : metrics) {
		n += metric.n;
		for (unsigned int i = 0u; i < metric.n; ++i) {
			if (aux::fcmp(min, metric.data[i], 0.00000001) > 0) {
				min = metric.data[i];
			}
			if (aux::fcmp(max, metric.data[i], 0.00000001) < 0) {
				max = metric.data[i];
			}
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

	for (auto&& metric : metrics) {
		for (unsigned int i = 0u; i < metric.n; ++i) {
			histogram.add(metric.data[i]);
		}
	}

	writer.write("{");
	writer.write("\"metric\":\"%s\",", metrics.front().name.c_str());
	writer.write("\"description\":\"%s\",", metrics.front().description.c_str());
	writer.write("\"n\":%d,", n);
	writer.write("\"min\":%g,\"max\":%g,", min, max);
	writer.write("\"hist\":[");
	vector<pair<double, unsigned int>> bins = histogram.get_bins();
	vector<pair<double, unsigned int>>::iterator bins_it = bins.begin();
	writer.write("%d", bins_it->second);
	if (bins.size() > 1) {
		while (++bins_it != bins.end() - 1) {
			writer.write(",%d", bins_it->second);
		}
	}
	writer.write("],");

	writer.write("\"mids\":[");
	if (bins.size() > 1) {
		bins_it = bins.begin() + 1;
		writer.write("%g", (bins_it->first - (bins_it->first - (bins_it - 1)->first) / 2.0));
		while (++bins_it < bins.end()) {
			writer.write(",%g", (bins_it->first - (bins_it->first - (bins_it - 1)->first) / 2.0));
		}
	} else {
		writer.write("%g", bins.front());
	}
	writer.write("]");
	writer.write("}");
}

