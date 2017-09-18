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
//#include "Histogram.h"
#include "TypeSwitcher.h"
#include "GzipWriter.h"
#include "aux.h"

#include <boost/program_options.hpp>

#include "Metric.h"
namespace po = boost::program_options;

using namespace std;
using namespace aux;

int main(int argc, char* argv[]) {
	vector<string> input_files;
	string metric_name("");
	unsigned int n_threads = 0u;
	string output_file("");

	po::options_description desc("Options");

	desc.add_options()
		("help,h", "Produce help message")
		("in,i", po::value<vector<string>>(&input_files)->multitoken()->required(), "Input VCF/BCF file(s). Must be indexed using tabix.")
		("metric,m", po::value<string>(&metric_name)->required(), "Metric name in INFO field.")
		("threads,t", po::value<unsigned int>(&n_threads)->required(), "Number of threads.")
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

	cout << n_threads << endl;

	try {
		vector<Metric> metrics;
		vector<thread> threads;

		for (auto&& input_file : input_files) {
			metrics.emplace_back(metric_name);
		}

		std::mutex cout_mutex;

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
//					cout << "Reading " << input_file << endl;

					cout_mutex.lock();
					cout << "Reading " << input_file << endl;
					cout_mutex.unlock();

//					metrics.emplace_back(metric_name);

					bcf_srs_t *sr = bcf_sr_init();

//					if (bcf_sr_add_reader(sr, input_file.c_str()) <= 0) {
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

					while (bcf_sr_next_line(sr) > 0) {
						bcf1_t* rec = bcf_sr_get_line(sr, 0);

						if ((rec->unpacked & BCF_UN_INFO) == 0) {
							bcf_unpack(rec, BCF_UN_INFO);
						}

			//			if ((rec->unpacked & BCF_UN_FLT) == 0) {
			//				bcf_unpack(rec, BCF_UN_FLT);
			//			}

			//			cout << rec->pos + 1 << " " << rec->d.id << " " << rec->n_info << " " << metrics_id << endl;

						for (int i = 0; i < rec->n_info; ++i) {
							if (rec->d.info[i].key == metrics_id) {
								uint8_t* info_ptr = rec->d.info[i].vptr;
								switch (rec->d.info[i].type) {
									case BCF_BT_INT8:
										for (int j = 0; j < rec->d.info[i].len; ++j) {
											int8_t value = le_to_i8(info_ptr + j * sizeof(int8_t));
											if (value != bcf_int8_missing) {
												metrics[f].add((double)value);
											}
										}
										break;
									case BCF_BT_INT16:
										for (int j = 0; j < rec->d.info[i].len; ++j) {
											int16_t value = le_to_i16(info_ptr + j * sizeof(int16_t));
											if (value != bcf_int16_missing) {
												metrics[f].add((double)value);
											}
										}
										break;
									case BCF_BT_INT32:
										for (int j = 0; j < rec->d.info[i].len; ++j) {
											int32_t value = le_to_i32(info_ptr + j * sizeof(int32_t));
											if (value != bcf_int32_missing) {
												metrics[f].add((double)value);
											}
										}
										break;
									case BCF_BT_FLOAT:
										for (int j = 0; j < rec->d.info[i].len; ++j) {
											uint32_t value = le_to_u32(info_ptr + j * sizeof(uint32_t));
											if (value != bcf_float_missing) {
												metrics[f].add((double)le_to_float(info_ptr + j * sizeof(uint32_t)));
											}
										}
										break;
									default:
										throw runtime_error("Error while resolving data types in INFO field!");
								}
							}
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


//		for (auto&& input_file : input_files) {
//			cout << "Reading " << input_file << endl;
//
//			metrics.emplace_back(metric_name);
//
//			bcf_srs_t *sr = bcf_sr_init();
//
//			if (bcf_sr_add_reader(sr, input_file.c_str()) <= 0) {
//				throw runtime_error("Error while initializing VCF/BCF reader!");
//			}
//
//			bcf_hdr_t* header = bcf_sr_get_header(sr, 0);
//
//			int metrics_id = bcf_hdr_id2int(header, BCF_DT_ID, metric_name.c_str());
//
//			if (bcf_hdr_idinfo_exists(header, BCF_HL_INFO, metrics_id) == 0) {
//				throw runtime_error("Specified INFO field was not found!");
//			}
//
//			int metrics_type = bcf_hdr_id2type(header, BCF_HL_INFO, metrics_id);
//			if ((metrics_type != BCF_HT_REAL) && (metrics_type != BCF_HT_INT)) {
//				throw runtime_error("Specified INFO field is not numeric!");
//			}
//
//			while (bcf_sr_next_line(sr) > 0) {
//				bcf1_t* rec = bcf_sr_get_line(sr, 0);
//
//				if ((rec->unpacked & BCF_UN_INFO) == 0) {
//					bcf_unpack(rec, BCF_UN_INFO);
//				}
//
//	//			if ((rec->unpacked & BCF_UN_FLT) == 0) {
//	//				bcf_unpack(rec, BCF_UN_FLT);
//	//			}
//
//	//			cout << rec->pos + 1 << " " << rec->d.id << " " << rec->n_info << " " << metrics_id << endl;
//
//				for (int i = 0; i < rec->n_info; ++i) {
//					if (rec->d.info[i].key == metrics_id) {
//						uint8_t* info_ptr = rec->d.info[i].vptr;
//						switch (rec->d.info[i].type) {
//							case BCF_BT_INT8:
//								for (int j = 0; j < rec->d.info[i].len; ++j) {
//									int8_t value = le_to_i8(info_ptr + j * sizeof(int8_t));
//									if (value != bcf_int8_missing) {
//										metrics.back().add((double)value);
//									}
//								}
//								break;
//							case BCF_BT_INT16:
//								for (int j = 0; j < rec->d.info[i].len; ++j) {
//									int16_t value = le_to_i16(info_ptr + j * sizeof(int16_t));
//									if (value != bcf_int16_missing) {
//										metrics.back().add((double)value);
//									}
//								}
//								break;
//							case BCF_BT_INT32:
//								for (int j = 0; j < rec->d.info[i].len; ++j) {
//									int32_t value = le_to_i32(info_ptr + j * sizeof(int32_t));
//									if (value != bcf_int32_missing) {
//										metrics.back().add((double)value);
//									}
//								}
//								break;
//							case BCF_BT_FLOAT:
//								for (int j = 0; j < rec->d.info[i].len; ++j) {
//									uint32_t value = le_to_u32(info_ptr + j * sizeof(uint32_t));
//									if (value != bcf_float_missing) {
//										metrics.back().add((double)le_to_float(info_ptr + j * sizeof(uint32_t)));
//									}
//								}
//								break;
//							default:
//								throw runtime_error("Error while resolving data types in INFO field!");
//						}
//					}
//				}
//			}
//
//			bcf_sr_destroy(sr);
//		}

		cout << "Computing histograms... " << flush;

		GzipWriter writer;

		writer.open(output_file.c_str());

//		metrics.write_histogram(writer, 40);
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
