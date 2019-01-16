#include <iostream>
#include <fstream>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <vector>
#include <string>
#include <set>
#include <thread>
#include <mutex>
#include "hts.h"
#include "vcf.h"
#include "synced_bcf_reader.h"
#include "TypeSwitcher.h"
#include <boost/program_options.hpp>
#include "Percentiles.h"
#include "aux.h"

namespace po = boost::program_options;

using namespace std;
using namespace aux;

inline void read_info_values(bcf_info_t &info, vector<double> &values) throw(runtime_error) {
    uint8_t *info_ptr = info.vptr;
    switch (info.type) {
        case BCF_BT_INT8:
            for (int j = 0; j < info.len; ++j) {
                int8_t value = le_to_i8(info_ptr + j * sizeof(int8_t));
                if (value != bcf_int8_missing) {
                    values.push_back((double) value);
                }
            }
            break;
        case BCF_BT_INT16:
            for (int j = 0; j < info.len; ++j) {
                int16_t value = le_to_i16(info_ptr + j * sizeof(int16_t));
                if (value != bcf_int16_missing) {
                    values.push_back((double) value);
                }
            }
            break;
        case BCF_BT_INT32:
            for (int j = 0; j < info.len; ++j) {
                int32_t value = le_to_i32(info_ptr + j * sizeof(int32_t));
                if (value != bcf_int32_missing) {
                    values.push_back((double) value);
                }
            }
            break;
        case BCF_BT_FLOAT:
            for (int j = 0; j < info.len; ++j) {
                uint32_t value = le_to_u32(info_ptr + j * sizeof(uint32_t));
                if (value != bcf_float_missing) {
                    values.push_back((double) le_to_float(info_ptr + j * sizeof(uint32_t)));
                }
            }
            break;
        default:
            throw runtime_error("Error while resolving data types in INFO field!");
    }
}

struct bcf_header_summary {
    int metrics_id;
    int af_id;
    int ac_id;
};

struct bcf_record_summary {
    vector<char *> alts;
    vector<double> ac_values;
    vector<double> af_values;
    vector<double> metric_values;
    float qual;
    unsigned int position;
    const char *chrom;
    char *id;
    char *ref;
    bool pass;

    bcf_record_summary() : pass(true) {}

    void clear() {
        alts.clear();
        ac_values.clear();
        af_values.clear();
        metric_values.clear();
        pass = true;
    }
};

void parse_bcf_header(bcf_hdr_t *header, string &metric_name, bcf_header_summary &h, bool qual) throw(runtime_error) {
    if (!qual) {
        h.metrics_id = bcf_hdr_id2int(header, BCF_DT_ID, metric_name.c_str());

        if (bcf_hdr_idinfo_exists(header, BCF_HL_INFO, h.metrics_id) == 0) {
            throw runtime_error("Specified INFO field was not found!");
        }

        int metrics_type = bcf_hdr_id2type(header, BCF_HL_INFO, h.metrics_id);
        if ((metrics_type != BCF_HT_REAL) && (metrics_type != BCF_HT_INT)) {
            throw runtime_error("Specified INFO field is not numeric!");
        }
    }

    h.af_id = bcf_hdr_id2int(header, BCF_DT_ID, "AF");
    if (bcf_hdr_idinfo_exists(header, BCF_HL_INFO, h.af_id) == 0) {
        throw runtime_error("The AF INFO field was not found!");
    }

    h.ac_id = bcf_hdr_id2int(header, BCF_DT_ID, "AC");
    if (bcf_hdr_idinfo_exists(header, BCF_HL_INFO, h.ac_id) == 0) {
        throw runtime_error("The AC INFO field was not found!");
    }
}

void parse_bcf_record(bcf_hdr_t *header, bcf1_t *rec, bcf_header_summary &h, bcf_record_summary &r,
                      bool qual) throw(runtime_error) {
    r.clear();
    if ((rec->unpacked & BCF_UN_INFO) == 0) {
        bcf_unpack(rec, BCF_UN_INFO);
    }
    if ((rec->unpacked & BCF_UN_FLT) == 0) {
        bcf_unpack(rec, BCF_UN_FLT);
    }
    for (int i = 0; i < rec->d.n_flt; ++i) {
        if (strcmp("PASS", bcf_hdr_int2id(header, BCF_DT_ID, rec->d.flt[i])) != 0) {
            r.pass = false;
            break;
        }
    }
    r.chrom = bcf_seqname(header, rec);
    r.position = rec->pos + 1;
    r.id = rec->d.id;
    r.ref = rec->d.allele[0];
    for (int i = 1; i < rec->n_allele; ++i) {
        r.alts.push_back(rec->d.allele[i]);
    }
    r.qual = rec->qual;

    for (int i = 0; i < rec->n_info; ++i) {
        if (rec->d.info[i].key == h.af_id) {
            read_info_values(rec->d.info[i], r.af_values);
        } else if (rec->d.info[i].key == h.ac_id) {
            read_info_values(rec->d.info[i], r.ac_values);
        } else if ((!qual) && (rec->d.info[i].key == h.metrics_id)) {
            read_info_values(rec->d.info[i], r.metric_values);
        }
    }

    if (r.ac_values.size() == 0 && r.af_values.size() == 0) {
        cerr << "[warning] No AC and AF INFO fields for " << r.chrom << ":" << r.position << ":" << r.ref << "/" << r.alts.at(0);
        for (unsigned int i = 1u; i < r.alts.size(); ++i) {
            cerr << "/" << r.alts.at(i);
        }
        cerr << ". Assuming AC = 0 and AF = 0.0 for all alternate alleles." << endl;
        for (unsigned int i = 0; i < r.alts.size(); ++i) {
            r.ac_values.push_back(0.0);
            r.af_values.push_back(0.0);
        }
    } else if ((r.alts.size() != r.ac_values.size()) || (r.ac_values.size() != r.af_values.size())) {
        
        cerr << "[warning] Number of values in AC and AF INFO fields don't match for " << r.chrom << ":" << r.position << ":" << r.ref << "/" << r.alts.at(0);
        for (unsigned int i = 1u; i < r.alts.size(); ++i) {
            cerr << "/" << r.alts.at(i);
        }
        cerr << ". Assuming AC = 0 and AF = 0.0 for all alternate alleles." << endl;
        for (unsigned int i = 0; i < r.alts.size(); ++i) {
            r.ac_values.push_back(0.0);
            r.af_values.push_back(0.0);
        }
    }
}


int main(int argc, char *argv[]) {
    vector<string> input_files;
    string metric_name("");
    unsigned int n_threads = 0u;
    string output_prefix("");
    double min_maf = -numeric_limits<double>::max();
    double max_maf = numeric_limits<double>::max();
    double allele_count = numeric_limits<double>::max();
    unsigned int n_percentiles = 0u;
    string metric_desc("");
    bool qual = false;

    po::options_description desc("Options");

    desc.add_options()
            ("help,h", "Produce help message")
            ("in,i", po::value<vector<string>>(&input_files)->multitoken()->required(),
             "Input VCF/BCF file(s). Must be indexed using tabix.")
            ("metric,m", po::value<string>(&metric_name), "Metric name in INFO field. By default, QUAL is used.")
            ("threads,t", po::value<unsigned int>(&n_threads)->required(), "Number of threads.")
            ("min-maf,f", po::value<double>(&min_maf), "Minimal minor allele frequency.")
            ("max-maf,F", po::value<double>(&max_maf), "Maximal minor allele frequency.")
            ("ac,a", po::value<double>(&allele_count), "Allele count (e.g. set to 1 to keep singletons only.")
            ("percentiles,p", po::value<unsigned int>(&n_percentiles)->default_value(100u), "Number of percentiles to compute. Default is 100.")
            ("desc,d", po::value<string>(&metric_desc), "Description (saved into JSON result).")
            ("out,o", po::value<string>(&output_prefix)->required(),
             "Prefix for output files. Output is compressed using bgzip.");

    po::variables_map vm;

    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        if (vm.count("help")) {
            cout
                    << "This program computes percentiles for the QUAL field or any arbitrary numeric INFO field across all variants."
                    << endl;
            cout << "Multi-allelic variants are treated as multiple bi-allelic variants." << endl << endl;
            cout << desc << endl;
            return 0;
        }
        po::notify(vm);
    } catch (po::error &e) {
        cout << "Error in command line:" << endl;
        cout << e.what() << endl;
        return 1;
    }

    vector<double> probabilities;
    for (double p = 0; p <= n_percentiles; p += 1) {
        probabilities.push_back(p / (double) n_percentiles);
    }

    if (n_threads > std::thread::hardware_concurrency()) {
        n_threads = std::thread::hardware_concurrency();
    }

    if (metric_name.length() == 0) {
        qual = true;
    }

    try {
        vector<Percentiles> percentiles_per_file;
        vector<thread> threads;

        std::mutex cout_mutex;

        for (auto it = input_files.begin(); it != input_files.end(); it++) {
            percentiles_per_file.emplace_back();
        }

        unsigned int step =
                n_threads > input_files.size() ? 1 : (unsigned int) lround(input_files.size() / (double) n_threads);
        unsigned int start = 0u, end = 0u;

        cout << "Reading files... " << endl;

        /* First pass to compute percentiles */
        while (end < input_files.size()) {
            if ((start + step) / step < n_threads) {
                end = start + step;
            } else {
                end = input_files.size();
            }

            threads.emplace_back([&](unsigned int start, unsigned int end) {
                for (unsigned int f = start; f < end; ++f) {
                    const char *input_file = input_files.at(f).c_str();

                    cout_mutex.lock();
                    cout << " " << input_file << endl;
                    cout_mutex.unlock();

                    bcf_header_summary header_summary;
                    bcf_record_summary record_summary;
                    double maf = 0.0;

                    bcf_srs_t *sr = bcf_sr_init();
                    if (bcf_sr_add_reader(sr, input_file) <= 0) {
                        throw runtime_error("Error while initializing VCF/BCF reader!");
                    }
                    bcf_hdr_t *header = bcf_sr_get_header(sr, 0);
                    parse_bcf_header(header, metric_name, header_summary, qual);
                    while (bcf_sr_next_line(sr) > 0) {
                        bcf1_t *rec = bcf_sr_get_line(sr, 0);
                        parse_bcf_record(header, rec, header_summary, record_summary, qual);
                        if (qual) {
                            for (unsigned int i = 0u; i < record_summary.alts.size(); ++i) {
                                if ((allele_count != numeric_limits<double>::max()) &&
                                    (allele_count != record_summary.ac_values[i])) {
                                    continue;
                                }
                                maf = fcmp(record_summary.af_values[i], 0.5, 0.00000001) > 0 ? 1.0 -
                                                                                               record_summary.af_values[i]
                                                                                             : record_summary.af_values[i];
                                if ((fcmp(maf, min_maf, 0.00000001) < 0) || (fcmp(maf, max_maf, 0.00000001) > 0)) {
                                    continue;
                                }
                                percentiles_per_file[f].add(record_summary.qual, record_summary.pass);
                            }
                        } else if (record_summary.metric_values.size() ==
                                   1) { // single metric value for all alt alleles
                            for (unsigned int i = 0u; i < record_summary.alts.size(); ++i) {
                                if ((allele_count != numeric_limits<double>::max()) &&
                                    (allele_count != record_summary.ac_values[i])) {
                                    continue;
                                }
                                maf = fcmp(record_summary.af_values[i], 0.5, 0.00000001) > 0 ? 1.0 -
                                                                                               record_summary.af_values[i]
                                                                                             : record_summary.af_values[i];
                                if ((fcmp(maf, min_maf, 0.00000001) < 0) || (fcmp(maf, max_maf, 0.00000001) > 0)) {
                                    continue;
                                }
                                percentiles_per_file[f].add(record_summary.metric_values[i], record_summary.pass);
                            }
                        } else if (record_summary.metric_values.size() ==
                                   record_summary.alts.size()) { // one metric value for every alt allele
                            for (unsigned int i = 0u; i < record_summary.alts.size(); ++i) {
                                if ((allele_count != numeric_limits<double>::max()) &&
                                    (allele_count != record_summary.ac_values[i])) {
                                    continue;
                                }
                                maf = fcmp(record_summary.af_values[i], 0.5, 0.00000001) > 0 ? 1.0 -
                                                                                               record_summary.af_values[i]
                                                                                             : record_summary.af_values[i];
                                if ((fcmp(maf, min_maf, 0.00000001) < 0) || (fcmp(maf, max_maf, 0.00000001) > 0)) {
                                    continue;
                                }
                                percentiles_per_file[f].add(record_summary.metric_values[i], record_summary.pass);
                            }
                        }
                    }
                    bcf_sr_destroy(sr);
                }
            }, start, end);

            start = end;
        }

        for (auto &&thread : threads) {
            thread.join();
        }

        cout << "Done." << endl;
        cout << "Merging... " << flush;
        Percentiles percentiles(qual ? "QUAL" : metric_name, qual ? "Phred-scaled site quality score" : metric_desc,
                                percentiles_per_file);
        cout << "Done." << endl;

        cout << "Computing percentiles... " << flush;
        percentiles.compute(probabilities);
        string output_percentiles_file = output_prefix + ".all_percentiles.json.gz";
        BGZF *ofp = bgzf_open(output_percentiles_file.c_str(), "w");
        if (!ofp) {
            throw runtime_error("Error while opening output file!");
        }
        percentiles.write(ofp);
        write(ofp, "\n");
        if (bgzf_close(ofp) != 0) {
            throw runtime_error("Error while closing output file!");
        }
        cout << "Done." << endl;

 
        /* Second pass to write percentiles */
        cout << "Writing variant percentiles... " << flush;
        string output_file = output_prefix + ".variant_percentile.vcf.gz";            
        BGZF *variant_ofp = bgzf_open(output_file.c_str(), "w");
        if (!variant_ofp) {
            throw runtime_error("Error while opening output file!");
        }
        write(variant_ofp, "##fileformat=VCFv4.3\n");
        if (qual) {
            write(variant_ofp, "##INFO=<ID=QUAL_PCTL,Number=2,Type=Float,Description=\"Phred-scaled quality score for the assertion made in ALT\">\n");
        } else {
            write(variant_ofp, "##INFO=<ID=%s_PCTL,Number=2,Type=Float,Description=\"%s\">\n", metric_name.c_str(), metric_desc.c_str());
        }
        write(variant_ofp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for (auto&& input_file : input_files) {
            bcf_header_summary header_summary;
            bcf_record_summary record_summary;
            double p_low = 0.0, p_high = 0.0;
            bcf_srs_t *sr = bcf_sr_init();
            if (bcf_sr_add_reader(sr, input_file.c_str()) <= 0) {
                throw runtime_error("Error while initializing VCF/BCF reader!");
            }
            bcf_hdr_t *header = bcf_sr_get_header(sr, 0);
            parse_bcf_header(header, metric_name, header_summary, qual);
            while (bcf_sr_next_line(sr) > 0) {
                bcf1_t *rec = bcf_sr_get_line(sr, 0);
                parse_bcf_record(header, rec, header_summary, record_summary, qual);
                if (qual) {
                    percentiles.probability(record_summary.qual, p_low, p_high);
                    for (unsigned int i = 0u; i < record_summary.alts.size(); ++i) {
                         write(variant_ofp, "%s\t%lu\t%s\t%s\t%s\t.\t.\tQUAL_PCTL=%g,%g\n", 
                               record_summary.chrom, record_summary.position, record_summary.id,
                               record_summary.ref, record_summary.alts[i], p_low, p_high);
                    }
                } else {
                    if (record_summary.metric_values.size() == 0) {
                         continue;
                    }
                    if (record_summary.metric_values.size() == 1) { // single metric for all alt alleles
                         percentiles.probability(record_summary.metric_values[0], p_low, p_high);
                         for (unsigned int i = 0u; i < record_summary.alts.size(); ++i) {
                              write(variant_ofp, "%s\t%lu\t%s\t%s\t%s\t.\t.\t%s_PCTL=%g,%g\n", 
                                    record_summary.chrom, record_summary.position, record_summary.id,
                                    record_summary.ref, record_summary.alts[i], metric_name.c_str(), p_low, p_high);
                         }
                    } else if (record_summary.metric_values.size() == record_summary.alts.size()) { // one metric for per alt allele
                         for (unsigned int i = 0u; i < record_summary.metric_values.size(); ++i) {
                             percentiles.probability(record_summary.metric_values[i], p_low, p_high);
                             write(variant_ofp, "%s\t%lu\t%s\t%s\t%s\t.\t.\t%s_PCTL=%g,%g\n", 
                                   record_summary.chrom, record_summary.position, record_summary.id,
                                   record_summary.ref, record_summary.alts[i], metric_name.c_str(), p_low, p_high);
                         }
                    } else {
                         throw runtime_error("Number of metrics doesn't match number of alternate alleles.");
                    }
                }
            }
            bcf_sr_destroy(sr);
        }
        if (bgzf_close(variant_ofp) != 0) {
            throw runtime_error("Error while closing output file!");
        }
        cout << "Done." << endl;
    } catch (exception &e) {
        cout << "Error: " << endl;
        cout << e.what() << endl;
        return 1;
    }

    return 0;
}
