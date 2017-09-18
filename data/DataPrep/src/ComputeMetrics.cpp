#include <iostream>
#include <fstream>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <vector>
#include <set>
#include "hts.h"
#include "vcf.h"
#include "synced_bcf_reader.h"
//#include "Histogram.h"
#include "TypeSwitcher.h"
#include "Metrics.h"
#include "GzipWriter.h"
#include "aux.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace std;
using namespace aux;

int main(int argc, char* argv[]) {
	vector<string> input_files;
	string metric_name("");
	string output_file("");

	po::options_description desc("Options");

	desc.add_options()
		("help,h", "Produce help message")
		("in,i", po::value<vector<string>>(&input_files)->multitoken()->required(), "Input VCF/BCF file(s). Must be indexed using tabix.")
		("metric,m", po::value<string>(&metric_name)->required(), "Metric name in INFO field.")
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

	try {
		Metrics metrics(metric_name);

		for (auto&& input_file : input_files) {
			cout << "Reading " << input_file << endl;

			bcf_srs_t *sr = bcf_sr_init();

			if (bcf_sr_add_reader(sr, input_file.c_str()) <= 0) {
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
										metrics.add((double)value);
									}
								}
								break;
							case BCF_BT_INT16:
								for (int j = 0; j < rec->d.info[i].len; ++j) {
									int16_t value = le_to_i16(info_ptr + j * sizeof(int16_t));
									if (value != bcf_int16_missing) {
										metrics.add((double)value);
									}
								}
								break;
							case BCF_BT_INT32:
								for (int j = 0; j < rec->d.info[i].len; ++j) {
									int32_t value = le_to_i32(info_ptr + j * sizeof(int32_t));
									if (value != bcf_int32_missing) {
										metrics.add((double)value);
									}
								}
								break;
							case BCF_BT_FLOAT:
								for (int j = 0; j < rec->d.info[i].len; ++j) {
									uint32_t value = le_to_u32(info_ptr + j * sizeof(uint32_t));
									if (value != bcf_float_missing) {
										metrics.add((double)le_to_float(info_ptr + j * sizeof(uint32_t)));
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

		cout << "Computing histograms... " << flush;

		GzipWriter writer;

		writer.open(output_file.c_str());

		metrics.write_histogram(writer, 40);

		writer.close();

		cout << "Done." << endl;
	} catch (exception &e) {
		cout << "Error: " << endl;
		cout << e.what() << endl;
      return 1;
   }

	return 0;
}
