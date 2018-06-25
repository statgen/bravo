#include "Histogram.h"

Histogram::Histogram(const vector<double>& borders) : total(0.0), n(0u), text(nullptr) {
	set<double> sorted_unique_borders;

	for (auto &&value : borders) {
		sorted_unique_borders.insert(value);
	}

	for (auto &&value : sorted_unique_borders) {
		this->bins.emplace_back(value, 0u);
	}

	last_added_value.first = 0.0;
	last_added_value.second = 0u;

	text = new char[1048576u];
}

Histogram::Histogram(Histogram &&histogram) {
	this->bins = std::move(histogram.bins);
	this->borders = std::move(histogram.borders);
	this->total = histogram.total;
	this->n = histogram.n;
	this->text = histogram.text;
	histogram.text = nullptr;
}

Histogram::~Histogram() noexcept {
	if (text != nullptr) {
		delete[] text;
		text = nullptr;
	}
}

void Histogram::add(double value) noexcept {
	total += value;
	++n;
	if (this->bins.size() == 1u) {
		if (aux::fcmp(value, this->bins[0].first, 0.00000001) == 0) {
			this->bins[0].second += 1u;
		}
	} else {
		last_added_value.first = value;
		auto it = upper_bound(this->bins.begin(), this->bins.end(), last_added_value, [] (const pair<double, unsigned int>& f, const pair<double, unsigned int>& s) -> bool {
			return (aux::fcmp(f.first, s.first, 0.00000001) < 0);
		});
		if (it != this->bins.begin()) {
			if (it != this->bins.end()) {
				(--it)->second +=1u;
			} else {
				if (aux::fcmp(value, (--it)->first,  0.00000001) <= 0) {
					(--it)->second += 1u;
				}
			}
		}
	}
}

vector<pair<double, unsigned int>> Histogram::get_bins() noexcept {
	vector<pair<double, unsigned int>> bins_copy;
	for (auto &&pair : bins) {
		bins_copy.emplace_back(pair.first, pair.second);
	}
	return bins_copy;
}

void Histogram::clear() noexcept {
	for (auto &&pair : bins) {
		pair.second = 0u;
	}
	total = 0.0;
	n = 0u;
}

double Histogram::get_total() noexcept {
	return total;
}

unsigned int Histogram::get_n() noexcept {
	return n;
}

double Histogram::get_average() noexcept {
	return total / (double)n;
}

const char* Histogram::get_text() throw (runtime_error) {
	int n = 0, i = 0;
	unsigned int bin = 1u;

	if ((n = sprintf(text, "%d", bins[0].second)) < 0) {
		throw runtime_error("Error while writing histogram!");
	}

	for (; bin < bins.size() - 1u; ++bin) {
		i += n;
		if ((n = sprintf(&text[i], "|%d", bins[bin].second)) < 0) {
			throw runtime_error("Error while writing histogram!");
		}
	}

	return text;
}


