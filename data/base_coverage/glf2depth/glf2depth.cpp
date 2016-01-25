#include <iostream>
#include <string>
#include <cstring>
#include <cstdint>
#include <fstream>
#include <memory>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <iomanip>
#include <exception>

#include <zlib.h>

using namespace std;

struct glf_cursor {
	uint32_t coordinate_offset;
	z_off_t glf_offset;

	glf_cursor(uint32_t coordinate_offset, z_off_t glf_offset) :  coordinate_offset(coordinate_offset), glf_offset(glf_offset) {

	}
};

struct stats {
	vector<unsigned int> depths;

	stats(unsigned int reserved_size) {
		depths.reserve(reserved_size);
	}
};

const unsigned int breaks[] = {1u, 5u, 10u, 15u, 20u, 25u, 30u, 50u, 100u};
const unsigned int n_breaks = 9u;

unsigned int n_samples = 0u;
string chromosome("");

unordered_set<string> chromosomes;
unordered_map<string, glf_cursor> glf_list;
unordered_map<unsigned int, stats> positions_map;
unordered_map<unsigned int, stats>::iterator positions_map_it;

void insert_depth(unsigned int position, unsigned int depth) {
	positions_map_it = positions_map.find(position);

	if (positions_map_it == positions_map.end()) {
		positions_map_it = positions_map.emplace(position, stats(n_samples)).first;
	}

	positions_map_it->second.depths.push_back(depth);
}

void write_depth(const string& output_file_name) throw (std::runtime_error) {
	ofstream output_stream;
	set<unsigned int> ordered_positions;

	unsigned int size = 0u;
	unsigned long int sum = 0ul;
	double mean = 0.0;
	double median = 0.0;
	unsigned int counts[n_breaks];

	output_stream.open(output_file_name, ios::binary | ios::app);
	if (output_stream.fail()) {
		throw std::runtime_error("Error while opening output file.");
	}

	for (auto&& position: positions_map) {
		ordered_positions.insert(position.first);
	}

	for (auto&& position: ordered_positions) {
		positions_map_it = positions_map.find(position);

		std::sort(positions_map_it->second.depths.begin(), positions_map_it->second.depths.end());

		size = positions_map_it->second.depths.size();
		if (size % 2 == 0) {
        	median = (positions_map_it->second.depths[size / 2 - 1] + positions_map_it->second.depths[size / 2]) / 2.0;
		} else {
			median =  positions_map_it->second.depths[(size - 1) / 2];
		}

		sum = 0ul;
		for (unsigned int i = 0; i < n_breaks; ++i) {
			counts[i] = 0u;
		}

		for (auto&& depth: positions_map_it->second.depths) {
			for (unsigned int i = n_breaks - 1u; i >= 0u; --i) {
				if (depth >= breaks[i]) {
					counts[i] += 1u;
					break;
				}
			}
			sum += depth;
		}

		for (int i = n_breaks - 2; i >= 0; --i) {
			counts[i] += counts[i + 1u];
		}

		mean = (double)sum / (double)size;

		output_stream << chromosome << "\t" << position << "\t";
		output_stream << "{\"chrom\":\"" << chromosome << "\",\"start\":" << position << ",\"end\":" << position << ",\"mean\":" << mean << ",\"median\":" << median;
		for (unsigned int i = 0u; i < n_breaks; ++i) {
			output_stream << ",\"" << breaks[i] << "\":" << ((double)counts[i] / (double)n_samples);
		}
		output_stream << "}" << endl;
	}

	if (output_stream.fail()) {
		throw std::runtime_error("Error while writing output file.");
	}

	output_stream.clear();
	output_stream.close();
	if (output_stream.fail()) {
		throw std::runtime_error("Error while closing output file.");
	}
}

void read_glf_list(const string& input_file_name, unordered_map<string, glf_cursor>& glf_list) throw (std::runtime_error) {
	ifstream glf_list_stream;
	unique_ptr<char[]> buffer = nullptr;
	unsigned int buffer_size = 160000u;

	glf_list_stream.open(input_file_name, ios::binary);
	if (glf_list_stream.fail()) {
		throw std::runtime_error("Error while opening GLF list.");
	}

	buffer = unique_ptr<char[]>(new char[buffer_size + 1u]);
	buffer[buffer_size] = '\0';

	while (!glf_list_stream.getline(buffer.get(), buffer_size).eof()) {
		glf_list.emplace(buffer.get(), glf_cursor(0, 0));
	}

	if (!glf_list_stream.eof() && glf_list_stream.fail()) {
		throw std::runtime_error("Error while reading GLF list.");
	}

	glf_list_stream.clear();
	glf_list_stream.close();
	if (glf_list_stream.fail()) {
		throw std::runtime_error("Error while closing GLF list.");
	}
}

void read_glf_header(const string& glf_name, string& sequence_name, z_off_t& glf_offset) throw (std::runtime_error) {
	gzFile gzfile = nullptr;
	int gzerrno = 0;

	int n_read_bytes = 0;

	char magic[4] = {'\0', '\0', '\0', '\0'};
	int32_t l_text = 0;
	int32_t l_name = 0;
	uint32_t ref_len = 0u;

	if (strcmp(zlibVersion(), ZLIB_VERSION) != 0) {
		throw std::runtime_error("Incompatible ZLIB versions.");
	}

	gzfile = gzopen(glf_name.c_str(), "rb");
	if (gzfile == nullptr) {
		throw std::runtime_error("Error while opening GLF.");
	}

//	Check magic field (char[4])
	if ((n_read_bytes = gzread(gzfile, magic, 4)) < 0) {
		throw std::runtime_error("Error while reading 'magic' field.");
	}

	if (magic[0] == 'G' && magic[1] == 'L' && magic[2] == 'F' && magic[3] == '\3') {
		cout << "- GLFv3 version detected." << endl;
	} else {
		throw std::runtime_error("File is not in GLFv3 format.");
	}

//	Check l_text field (int32_t)
	if (((n_read_bytes = gzread(gzfile, &l_text, sizeof(l_text))) < 0) || (n_read_bytes < (long int)sizeof(l_text))) {
		throw std::runtime_error("Error while reading 'l_text' field.");
	}
	cout << "- GLFv3 header length: " << l_text << endl;

//	Check text field (char[l_text])
	if (l_text > 0) {
		unique_ptr<char[]> text = unique_ptr<char[]>(new char[l_text]);
		if (((n_read_bytes = gzread(gzfile, text.get(), l_text)) < 0) || (n_read_bytes < l_text)) {
			throw std::runtime_error("Error while reading 'text' field.");
		}
		cout << "- GLFv3 header: " << text.get() << endl;
	}

//	Check l_name field (int32_t)
	n_read_bytes =  gzread(gzfile, &l_name, sizeof(l_name));
    if (n_read_bytes == 0) {
    	glf_offset = 0;
    } else {
    	if (n_read_bytes < (long int)sizeof(l_name)) {
    		throw std::runtime_error("Error while reading 'l_name' field.");
    	}
    	cout << "- Length of the reference sequence name: " << l_name << endl;

//		Check name field (char[l_name])
    	if (l_name > 0) {
    		unique_ptr<char[]> name = unique_ptr<char[]>(new char[l_name]);
    		if (((n_read_bytes = gzread(gzfile, name.get(), l_name)) < 0) || (n_read_bytes < l_name)) {
    			throw std::runtime_error("Error while reading 'name' field.");
    		}
    		cout << "- Reference sequence name: " << name.get() << endl;
    		sequence_name = name.get();
    	}

//		Check ref_len field (uint32_t)
    	if (((n_read_bytes = gzread(gzfile, &ref_len, sizeof(ref_len))) < 0) || (n_read_bytes < (long int)sizeof(ref_len))) {
    		throw std::runtime_error("Error while reading 'ref_len' field.");
    	}
    	cout << "- Length of the reference sequence: " << ref_len << endl;

    	glf_offset = gztell(gzfile);
    }

	gzerrno = gzclose(gzfile);
	gzfile = nullptr;
	if (gzerrno != Z_OK) {
		throw std::runtime_error("Error while closing GLF.");
	}
}

bool read_glf_entries(const string& glf_name, unsigned int max_position, uint32_t& coordinate_offset, z_off_t& glf_offset) throw (std::runtime_error) {
	gzFile gzfile = nullptr;
	int gzerrno = 0;

	int n_read_bytes = 0;
	bool eof = false;

	uint8_t rtype_ref = 0u;
	int rtype = 0;
	uint32_t offset = 0u;
	uint32_t min_depth = 0u;
	unsigned int read_depth = 0u;
	unsigned int coordinate = 0u;
	int16_t indelLen1 = 0;
	int16_t indelLen2 = 0;

	gzfile = gzopen(glf_name.c_str(), "rb");
	if (gzfile == nullptr) {
		throw std::runtime_error("Error while opening GLF.");
	}

	if (gzseek(gzfile, glf_offset, SEEK_SET) < 0) {
		throw std::runtime_error("Error while seeking the offset.");
	}

	coordinate = coordinate_offset;

	while (true) {
		glf_offset = gztell(gzfile);
		coordinate_offset = coordinate;

//		Check rtype_ref (uint8_t)
		if (((n_read_bytes = gzread(gzfile, &rtype_ref, sizeof(rtype_ref))) == 0)) {
			eof = true;
			break;
		}

		if ((n_read_bytes < 0) || (n_read_bytes < (long int)sizeof(rtype_ref))) {
			throw std::runtime_error("Error while reading 'rtype_ref' field.");
		}

		rtype = rtype_ref >> 4;

		if (rtype == 0) {
			eof = true;
			break;
		} else if ((rtype != 1) && (rtype != 2)) {
			throw std::runtime_error("Unrecognized record type.");
		} else {
//			Check offset (uint32_t)
			if (((n_read_bytes = gzread(gzfile, &offset, sizeof(offset))) < 0) || (n_read_bytes < (long int)sizeof(offset))) {
				throw std::runtime_error("Error while reading 'offset' field.");
			}

//			Check min_depth (uint32_t)
			if (((n_read_bytes = gzread(gzfile, &min_depth, sizeof(min_depth))) < 0) || (n_read_bytes < (long int)sizeof(min_depth))) {
				throw std::runtime_error("Error while reading 'min_depth' field.");
			}

			read_depth = min_depth & 0xFFFFFF;
			coordinate = coordinate + offset;

			if (coordinate > max_position) {
				eof = false;
				break;
			}

			if (offset != 0) {
				insert_depth(coordinate + 1u, read_depth);
			}

			if (rtype == 1) {
//				Skip rmsMapQ (uint8_t) and lk fields (uint8_t[10])
				if (gzseek(gzfile, 11 * sizeof(uint8_t), SEEK_CUR) < 0) {
					throw std::runtime_error("Error while seeking the end of the simple record.");
				}
			} else if (rtype == 2) {
//				Skip rmsMapQ (uint8_t), lkHom1 (uint8_t), lkHom2 (uint8_t), and lkHet (uint8_t) fields;
				if (gzseek(gzfile, 4 * sizeof(uint8_t), SEEK_CUR) < 0) {
					throw std::runtime_error("Error while seeking inside the indel record.");
				}

//				Check indelLen1 (int16_t)
				if (((n_read_bytes = gzread(gzfile, &indelLen1, sizeof(indelLen1))) < 0) || (n_read_bytes < (long int)sizeof(indelLen1))) {
					throw std::runtime_error("Error while reading 'indelLen1' field.");
				}

//				Check indelLen2 (int16_t)
				if (((n_read_bytes = gzread(gzfile, &indelLen2, sizeof(indelLen2))) < 0) || (n_read_bytes < (long int)sizeof(indelLen2))) {
					throw std::runtime_error("Error while reading 'indelLen2' field.");
				}

				if (indelLen1 < 0) {
					indelLen1 = abs(indelLen1);
				}

				if (indelLen2 < 0) {
					indelLen2 = abs(indelLen2);
				}

//				Skip indelSeq1 (char[indelLen1]) and indelSeq2 (char[indelLen2]) fields
				if (gzseek(gzfile, (indelLen1 + indelLen2) * sizeof(char), SEEK_CUR) < 0) {
					throw std::runtime_error("Error while seeking the end of the indel record.");
				}
			}
		}
	}

	gzerrno = gzclose(gzfile);
	gzfile = nullptr;
	if (gzerrno != Z_OK) {
		throw std::runtime_error("Error while closing GLF.");
	}

	return eof;
}

/*
 *
 * Description:
 * 	Constructs depth histograms for every base in provided GLF files.
 *
 * Arguments:
 * 	1 -- List of GLF files. One file per line.
 * 	2 -- Output VCF name.
 * 	3 -- Chunk size in base pairs e.g. 1000000.
 *
 */
int main(int argc, char* argv[]) {
	string glf_list_file(argv[1]);
	string output_depth_file(argv[2]);
	string sequence_name("");
	unsigned long int chunk = 0u;
	unsigned long int max_position = 0u;
	unsigned int iteration = 0u;
	bool eof = false;
	bool eof_all = false;

	try {
		chunk = std::stoul(argv[3]);

		read_glf_list(glf_list_file, glf_list);

		cout << "HEADERS" << endl;
		for (auto&& glf : glf_list) {
			cout << "Processing " << glf.first << endl;
			read_glf_header(glf.first, sequence_name, glf.second.glf_offset);
			if (glf.second.glf_offset != 0) {
				chromosomes.emplace(sequence_name);
				++n_samples;
			}
			cout << "Done" << endl;
		}
		cout << endl;

		if (chromosomes.size() > 1) {
			throw std::runtime_error("Multiple chromosomes are not permitted.");
		} else if (chromosomes.size() == 0) {
			throw std::runtime_error("All GLF files are empty.");
		}
		chromosome = *(chromosomes.begin());

		while (!eof_all) {
			eof_all = true;
			++iteration;
			max_position += chunk;
			cout << "ITERATION "<< iteration << " (<= " << max_position << " bp)" << endl;
			for (auto&& glf : glf_list) {
				if (glf.second.glf_offset > 0) {
					cout << "Processing " << glf.first << " (offset = <" << glf.second.coordinate_offset  << ", " << glf.second.glf_offset << ">)" << endl;
					eof = read_glf_entries(glf.first, max_position, glf.second.coordinate_offset, glf.second.glf_offset);
					cout << "Done (offset = <" << glf.second.coordinate_offset  << ", " << glf.second.glf_offset << ">)" << endl;
				}
			}

			cout << "Writing to file " << output_depth_file << endl;
			write_depth(output_depth_file);
			positions_map.clear();
			cout << "Done" << endl;

			eof_all = eof_all &&  eof;
		}
	} catch (std::exception& e) {
		cerr << e.what() << endl;
	}

	return 0;
}
