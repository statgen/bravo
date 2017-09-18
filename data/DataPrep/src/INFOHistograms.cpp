#include <iostream>
#include <fstream>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <vector>
#include <set>
#include <thread>
#include <mutex>
#include "hts.h"
#include "vcf.h"
#include "synced_bcf_reader.h"
#include "TypeSwitcher.h"
#include "GzipWriter.h"
#include "aux.h"

#include <boost/program_options.hpp>

#include "Metric.h"
namespace po = boost::program_options;

using namespace std;
using namespace aux;

inline void read_info_values(bcf_info_t &info, vector<double>& values) throw (runtime_error) {
	uint8_t* info_ptr = info.vptr;
	switch (info.type) {
	case BCF_BT_INT8:
		for (int j = 0; j < info.len; ++j) {
			int8_t value = le_to_i8(info_ptr + j * sizeof(int8_t));
			if (value != bcf_int8_missing) {
				values.push_back((double)value);
			}
		}
		break;
	case BCF_BT_INT16:
		for (int j = 0; j < info.len; ++j) {
			int16_t value = le_to_i16(info_ptr + j * sizeof(int16_t));
			if (value != bcf_int16_missing) {
				values.push_back((double)value);
			}
		}
		break;
	case BCF_BT_INT32:
		for (int j = 0; j < info.len; ++j) {
			int32_t value = le_to_i32(info_ptr + j * sizeof(int32_t));
			if (value != bcf_int32_missing) {
				values.push_back((double)value);
			}
		}
		break;
	case BCF_BT_FLOAT:
		for (int j = 0; j < info.len; ++j) {
			uint32_t value = le_to_u32(info_ptr + j * sizeof(uint32_t));
			if (value != bcf_float_missing) {
				values.push_back((double)le_to_float(info_ptr + j * sizeof(uint32_t)));
			}
		}
		break;
	default:
		throw runtime_error("Error while resolving data types in INFO field!");
	}
}

int main(int argc, char* argv[]) {
	vector<string> input_files;
	string metric_name("");
	unsigned int n_threads = 0u;
	string output_file("");
	double min_maf = -numeric_limits<double>::max();
	double max_maf = numeric_limits<double>::max();
	double allele_count = numeric_limits<double>::max();
	string metric_desc("");

	po::options_description desc("Options");

	desc.add_options()
		("help,h", "Produce help message")
		("in,i", po::value<vector<string>>(&input_files)->multitoken()->required(), "Input VCF/BCF file(s). Must be indexed using tabix.")
		("metric,m", po::value<string>(&metric_name)->required(), "Metric name in INFO field.")
		("threads,t", po::value<unsigned int>(&n_threads)->required(), "Number of threads.")
		("min-maf,f", po::value<double>(&min_maf), "Minimal minor allele frequency.")
		("max-maf,F", po::value<double>(&max_maf), "Maximal minor allele frequency.")
		("ac,a", po::value<double>(&allele_count), "Allele count (e.g. set to 1 to keep singletons only.")
		("desc,d", po::value<string>(&metric_desc), "Description (saved into JSON result).")
		("out,o", po::value<string>(&output_file)->required(), "Output file. Compressed using gzip.")
	;

	po::variables_map vm;

	try {
		po::store(po::parse_command_line(argc, argv, desc), vm);
		if (vm.count("help")) {
			cout << "This program computes histogram for the specified INFO field across all variants." << endl << endl;
			cout << desc << endl;
			return 0;
		}
		po::notify(vm);
	} catch (po::error &e) {
		cout << "Error in command line:" << endl;
		cout << e.what() << endl;
		return 1;
	}

	if (n_threads > std::thread::hardware_concurrency()) {
		n_threads = std::thread::hardware_concurrency();
	}

	try {
		vector<Metric> metrics;
		vector<thread> threads;
		std::mutex cout_mutex;

		for (auto&& input_file : input_files) {
			metrics.emplace_back(metric_name, metric_desc);
		}

		unsigned int step = n_threads > input_files.size() ? 1 : (unsigned int)lround(input_files.size() / (double)n_threads);
		unsigned int start = 0u, end = 0u;
		while (end < input_files.size()) {
			if ((start + step) / step < n_threads) {
				end = start + step;
			} else {
				end = input_files.size();
			}

			threads.emplace_back([&](unsigned int start, unsigned int end) {
				for (unsigned int f = start; f < end; ++f) {
					const char* input_file = input_files.at(f).c_str();

					cout_mutex.lock();
					cout << "Reading " << input_file << endl;
					cout_mutex.unlock();

					bcf_srs_t *sr = bcf_sr_init();

					if (bcf_sr_add_reader(sr, input_file) <= 0) {
						throw runtime_error("Error while initializing VCF/BCF reader!");
					}

					bcf_hdr_t* header = bcf_sr_get_header(sr, 0);

					int metrics_id = bcf_hdr_id2int(header, BCF_DT_ID, metric_name.c_str());

					if (bcf_hdr_idinfo_exists(header, BCF_HL_INFO, metrics_id) == 0) {
						throw runtime_error("Specified INFO field was not found!");
					}

					int metrics_type = bcf_hdr_id2type(header, BCF_HL_INFO, metrics_id);
					if ((metrics_type != BCF_HT_REAL) && (metrics_type != BCF_HT_INT)) {
						throw runtime_error("Specified INFO field is not numeric!");
					}

					int af_id = bcf_hdr_id2int(header, BCF_DT_ID, "AF");
					if (bcf_hdr_idinfo_exists(header, BCF_HL_INFO, af_id) == 0) {
						throw runtime_error("The AF INFO field was not found!");
					}

					int ac_id = bcf_hdr_id2int(header, BCF_DT_ID, "AC");
					if (bcf_hdr_idinfo_exists(header, BCF_HL_INFO, ac_id) == 0) {
						throw runtime_error("The AC INFO field was not found!");
					}

					vector<double> ac_values;
					vector<double> af_values;
					vector<double> metric_values;
					double maf = 0.0;

					while (bcf_sr_next_line(sr) > 0) {
						bcf1_t* rec = bcf_sr_get_line(sr, 0);

						if ((rec->unpacked & BCF_UN_INFO) == 0) {
							bcf_unpack(rec, BCF_UN_INFO);
						}

			//			if ((rec->unpacked & BCF_UN_FLT) == 0) {
			//				bcf_unpack(rec, BCF_UN_FLT);
			//			}

			//			cout << rec->pos + 1 << " " << rec->d.id << " " << rec->n_info << " " << metrics_id << endl;

						ac_values.clear();
						af_values.clear();
						metric_values.clear();

						for (int i = 0; i < rec->n_info; ++i) {
							bcf_info_t info = rec->d.info[i];

							if (rec->d.info[i].key == af_id) {
								read_info_values(rec->d.info[i], af_values);
							} else if (rec->d.info[i].key == ac_id) {
								read_info_values(rec->d.info[i], ac_values);
							} else if (rec->d.info[i].key == metrics_id) {
								read_info_values(rec->d.info[i], metric_values);
							}
						}

						if ((ac_values.size() != af_values.size()) && (ac_values.size() != metric_values.size())) {
							throw runtime_error("Missing values in AC, AF and the specified INFO fields!");
						}

						for (unsigned int i = 0u; i < metric_values.size(); ++i) {
							if ((allele_count != numeric_limits<double>::max()) && (allele_count != ac_values[i])) {
								continue;
							}
							maf = fcmp(af_values[i], 0.5, 0.00000001) > 0 ? 1.0 - af_values[i] : af_values[i];
							if ((fcmp(maf, min_maf, 0.00000001) < 0) || (fcmp(maf, max_maf, 0.00000001) > 0)) {
								continue;
							}
							metrics[f].add(metric_values[i]);
						}
					}

					bcf_sr_destroy(sr);
				}
			}, start, end);

			start = end;
		}

		for (auto&& thread : threads) {
			thread.join();
		}

		cout << "Computing histograms... " << flush;

		GzipWriter writer;
		writer.open(output_file.c_str());
		Metric::write_histogram(writer, metrics, 40);
		writer.close();

		cout << "Done." << endl;
	} catch (exception &e) {
		cout << "Error: " << endl;
		cout << e.what() << endl;
      return 1;
   }

	return 0;
}
