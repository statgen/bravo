#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <boost/program_options.hpp>
#include <vector>
#include <regex>
#include <unordered_map>
#include <tuple>
#include <unordered_set>
#include <chrono>
#include <thread>
#include "hts.h"
#include "synced_bcf_reader.h"
#include "htslib/sam.h"

namespace po = boost::program_options;

using namespace std;

struct Variant {
    unsigned int pos;
    string chrom;
    string ref;
    string alt;
    Variant(const string& chrom, unsigned int pos, const string& ref, const string& alt) : chrom(chrom), pos(pos), ref(ref), alt(alt) {}
};

struct Region {
    unsigned int variant_idx;
    unsigned int sample_idx;
    bool is_hom;
    Region(unsigned int variant_idx, unsigned int sample_idx, bool is_hom): variant_idx(variant_idx), sample_idx(sample_idx), is_hom(is_hom) {}
};

vector<string> optional_tag_names {
    "AM", "AS", "BC", "BQ", "BZ", "CB",
    "CC", "CG", "CM", "CO", "CP", "CQ",
    "CR", "CS", "CT", "CY", "E2", "FI",
    "FS", "FZ", "GC", "GQ", "GS", "H0",
    "H1", "H2", "HI", "IH", "LB", "MC",
    "MD", "MF", "MI", "MQ", "NH", "NM",
    "OA", "OC", "OP", "OQ", "OX", "PG",
    "PQ", "PT", "PU", "Q2", "QT", "QX",
    "R2", "RG", "RT", "RX", "S2", "SA",
    "SM", "SQ", "TC", "TS", "U2", "UQ",
    "XA", "XS"
};

samFile* open_cram_with_retry(const char* cram_path, unsigned int max_retries = 10, unsigned int sleep_seconds = 60) {
    unsigned int retry = 1u;
    samFile *in_cram_fp = sam_open(cram_path, "r");
    while ((!in_cram_fp) && (retry <= max_retries)) {
        cout << "Re-trying (" << retry << ") to open " << cram_path << endl << flush;
        this_thread::sleep_for(chrono::seconds(sleep_seconds));
        in_cram_fp = sam_open(cram_path, "r");
        ++retry;
    }
    return in_cram_fp;
}

hts_idx_t* open_index_with_retry(samFile* in_cram_fp, const char* cram_path, const char* crai_path, unsigned int max_retries = 10, unsigned int sleep_seconds = 60) {
    unsigned int retry = 1u;
    hts_idx_t *in_idx = sam_index_load2(in_cram_fp, cram_path, crai_path);
    while ((!in_idx) && (retry <= max_retries)) {
        cout << "Re-trying (" << retry << ") to open " << cram_path << endl << flush;
        this_thread::sleep_for(chrono::seconds(sleep_seconds));
        in_idx = sam_index_load2(in_cram_fp, cram_path, crai_path);
        ++retry;
    }
    return in_idx;
}

void process_sample(
        const string& cram_path, const string& crai_path, const string& reference_path,
        vector<Region>& regions, samFile* out_cram_fp, bam_hdr_t* out_header, vector<Variant>& variants, unsigned int& window) {
    samFile *in_cram_fp = open_cram_with_retry(cram_path.c_str());
    if (!in_cram_fp) {
        throw runtime_error("Error while opening input CRAM file.");
    }
    if (hts_set_fai_filename(in_cram_fp, reference_path.c_str()) != 0) {
        throw runtime_error("Error while opening reference FASTA file.");
    }
    bam_hdr_t *in_header = sam_hdr_read(in_cram_fp);
    if (!in_header) {
        throw runtime_error("Error while reading header from input CRAM file.");
    }
    hts_idx_t *in_idx = open_index_with_retry(in_cram_fp, cram_path.c_str(), crai_path.c_str());
    if (!in_idx) {
        throw runtime_error("Error while reading input CRAI file.");
    }
    unsigned int start = 0u, stop = 0u;
    int ret = 0;
    const char* qname = nullptr;
    unsigned int qname_idx  = 0u;
    uint8_t *tag = nullptr;
    string new_qname("");
    bam1_t *b = bam_init1();
    if (b == nullptr) {
        throw runtime_error("Error allocating BAM structure.");
    }
    for (auto &&region : regions) {
        Variant& variant = variants.at(region.variant_idx);
        int tid = bam_name2id(in_header, variant.chrom.c_str());
        if (tid < 0) {
            throw runtime_error("Error while getting reference id from input CRAM file.");
        }
        start = window > variant.pos ? 0u : variant.pos - window;
        stop = variant.pos + window;

        hts_itr_t *iter = sam_itr_queryi(in_idx, tid, start, stop); // zero-based
        if (iter == nullptr) {
            throw runtime_error("Error while creating CRAM iterator.");
        }
        unordered_map<string, unsigned int> qnames;
        auto qnames_it = qnames.begin();
        while ((ret = sam_itr_next(in_cram_fp, iter, b)) >= 0) {
            qname = bam_get_qname(b);
            qnames_it = qnames.find(qname);
            if (qnames_it == qnames.end()) {
                qname_idx = qnames.emplace(qname, qnames.size() + 1).first->second;
            } else {
                qname_idx = qnames_it->second;
            }
            new_qname = to_string(variant.pos) + ":" + variant.ref + ":" + variant.alt + ":" + (region.is_hom ? "0" : "") + to_string(region.sample_idx) + ":" + to_string(qname_idx);
            if (bam_set_qname(b, new_qname.c_str()) < 0) {
                throw runtime_error("Error while setting new QNAME.");
            }
            for (auto &&tag_name: optional_tag_names) {
                if ((tag = bam_aux_get(b, tag_name.c_str())) != nullptr) {
                    if (bam_aux_del(b, tag) < 0) {
                        throw runtime_error("Error while deleting optional fields from CRAM file.");
                    }
                }
            }
            if (sam_write1(out_cram_fp, out_header, b) < 0) {
                throw runtime_error("Error while writing sequence read to output CRAM file.");
            }
        }
        hts_itr_destroy(iter);
    }

    bam_destroy1(b);
    sam_hdr_destroy(in_header);
    hts_idx_destroy(in_idx);
    sam_close(in_cram_fp);
}

int main(int argc, char* argv[]) {
    string input_file("");
    string crams_file("");
    string reference_file("");
    string output_cram_file("");
    unsigned int window = 0u;

    po::options_description desc("Options");

    desc.add_options()
            ("help,h", "Produce help message")
            ("in,i", po::value<string>(&input_file)->required(), "Input compressed VCF/BCF with HET and HOM INFO fields. Multi-allelic variants must be split into bi-allelic entries.")
            ("crams,c", po::value<string>(&crams_file)->required(), "Input file with a sample name and the corresponding CRAM and CRAI file paths per line. No header. Tab separated.")
            ("reference,r", po::value<string>(&reference_file)->required(), "Reference FASTA file.")
            ("window,w", po::value<unsigned int>(&window)->default_value(100), "Window size around each variant in base-pairs. Default 100 bp.")
            ("out,o", po::value<string>(&output_cram_file)->required(), "Ouptut CRAM file.")
            ;

    po::variables_map vm;

    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        if (vm.count("help")) {
            cout << "This program generates CRAM file with sequences from heterozygous/homozygous samples. CRAM file is not sorted and is compressed without a reference (i.e. 'no_ref' option)."  << endl << endl;
            cout << desc << endl;
            return 0;
        }
        po::notify(vm);
    } catch (po::error &e) {
        cout << "Error in command line:" << endl;
        cout << e.what() << endl;
        return 1;
    } catch (invalid_argument &e) {
        cout << "Error in command line:" << endl;
        cout << e.what() << endl;
        return 1;
    }

    // BEGIN: Read list of samples with corresponding CRAM and CRAI files.
    unordered_map<string, tuple<string, string> > crams;
    ifstream if_crams_file(crams_file);
    string line;
    vector<string> tokens;
    auto separator = regex("[ \t]");
    while (getline(if_crams_file, line)) {
        copy(sregex_token_iterator(line.begin(), line.end(), separator, -1), sregex_token_iterator(), back_inserter(tokens));
        if (tokens.size() != 3) {
            throw runtime_error("File with CRAM and CRAI files must have 3 tab-delimited fields.");
        }
        crams.emplace(tokens.at(0), make_tuple(tokens.at(1), tokens.at(2)));
        tokens.clear();
    }
    if (crams.size() == 0) {
        cout << "No samples (CRAM and CRAI files) listed." << endl;
        return 0;
    }
    // END.

    // BEGIN: Read VCF with HET and HOM variant carriers
    unordered_set<string> chrom;
    unsigned int pos = 0u;
    unsigned int start = numeric_limits<unsigned int>::max();
    unsigned int stop = 0u;
    unsigned int max_hom = 0u;
    unsigned int max_het = 0u;
    vector<Variant> variants;
    unsigned int variant_idx = 0u;
    unordered_map<string, vector<Region>> samples;

    bcf_srs_t *sr = bcf_sr_init();
    if (bcf_sr_add_reader(sr, input_file.c_str()) <= 0) {
        throw runtime_error("Error while initializing VCF/BCF reader!");
    }
    bcf_hdr_t* header = bcf_sr_get_header(sr, 0);
    int info_het_id = bcf_hdr_id2int(header, BCF_DT_ID, "HET");
    if (bcf_hdr_idinfo_exists(header, BCF_HL_INFO, info_het_id) == 0) {
        throw runtime_error("The HET INFO field was not found!");
    }
    int info_hom_id = bcf_hdr_id2int(header, BCF_DT_ID, "HOM");
    if (bcf_hdr_idinfo_exists(header, BCF_HL_INFO, info_hom_id) == 0) {
        throw runtime_error("The HOM INFO field was not found!");
    }
    while (bcf_sr_next_line(sr) > 0) {
        bcf1_t *rec = bcf_sr_get_line(sr, 0);
        if ((rec->unpacked & BCF_UN_INFO) == 0) {
            bcf_unpack(rec, BCF_UN_INFO);
        }
        if (rec->n_allele != 2) {
            throw runtime_error("Variant with more than 1 alternate allele in a single VCF record was detected. Only bi-allelic records are allowed.");
        }
        chrom.emplace(bcf_seqname(header, rec));
        if (chrom.size() > 1) {
            throw runtime_error("More than 1 chromosome detected in VCF. Only single chromosome per file is allowed.");
        }
        pos = rec->pos + 1;
        if (pos < start) {
            start = pos;
        }
        if (pos > stop) {
            stop = pos;
        }
        variants.emplace_back(bcf_seqname(header, rec), pos, rec->d.allele[0], rec->d.allele[1]);
        variant_idx = variants.size() - 1;
        for (int i = 0; i < rec->n_info; ++i) {
            if (rec->d.info[i].key == info_het_id) {
                if (rec->d.info[i].type != BCF_BT_CHAR) {
                    throw runtime_error("The HET INFO field is not a string.");
                }
                string sample("");
                for (int j = 0, idx = 1; j < rec->d.info[i].len; ++j) {
                    if (rec->d.info[i].vptr[j] != ',') {
                        sample.push_back(rec->d.info[i].vptr[j]);
                        if (j < rec->d.info[i].len - 1) {
                            continue;
                        }
                    }
                    samples.emplace(sample, vector<Region>()).first->second.emplace_back(variant_idx, idx, false);
                    if (idx > max_het) {
                        max_het = idx;
                    }
                    sample.clear();
                    ++idx;
                }
            } else if (rec->d.info[i].key == info_hom_id) {
                if (rec->d.info[i].type != BCF_BT_CHAR) {
                    throw runtime_error("The HOM INFO field is not a string.");
                }
                string sample("");
                for (int j = 0, idx = 1; j < rec->d.info[i].len; ++j) {
                    if (rec->d.info[i].vptr[j] != ',') {
                        sample.push_back(rec->d.info[i].vptr[j]);
                        if (j < rec->d.info[i].len - 1) {
                            continue;
                        }
                    }
                    samples.emplace(sample, vector<Region>()).first->second.emplace_back(variant_idx, idx, true);
                    if (idx > max_hom) {
                        max_hom = idx;
                    }
                    sample.clear();
                    ++idx;
                }
            }
        }
    }
    bcf_sr_destroy(sr);
    if (variants.size() == 0) {
        cout << "No variants in input VCF/BCF." << endl;
        return 0;
    }
    // END.

    // BEGIN: Sanity check to see if all randomly selected samples have corresponding CRAM/CRAI files.
    for (auto &&sample: samples) {
        if (crams.find(sample.first) == crams.end()) {
            throw runtime_error("Sample(s) without CRAM and CRAI files detected.");
        }
    }
    // END.

    // BEGIN: Open output CRAM and create output header.
    htsFormat cram_fmt;
    memset(&cram_fmt, 0, sizeof(cram_fmt));
    if (hts_parse_format(&cram_fmt, "cram,no_ref") < 0) { // disable reference-based comression with no_ref since we write out non-sorted data and don't want heavy IO/memory load when compressing.
        throw runtime_error("Error while setting format parameters for output CRAM file.");
    }
    samFile *out_cram_fp = hts_open_format(output_cram_file.c_str(), "wc", &cram_fmt);
    if (!out_cram_fp) {
        throw runtime_error("Error while opening output CRAM file.");
    }
    sam_hdr_t *out_header = sam_hdr_init();
    if (out_header == nullptr) {
        throw runtime_error("Error while creating header for output CRAM file.");
    }
    if (sam_hdr_add_line(out_header, "HD", "VN", "1.6", "SO", "unsorted", nullptr) != 0) {
        throw runtime_error("Error while writing 'HD' header line to output CRAM file.");
    }
    string co_value("MAX_HOM=" + to_string(max_hom) + ";MAX_HET=" + to_string(max_het));
    if (sam_hdr_add_line(out_header, "CO", co_value.c_str(), nullptr) != 0) {
        throw runtime_error("Error while writing 'HD' header line to output CRAM file.");
    }
    co_value = "REGION=" + *begin(chrom)  + ":" + to_string(start) + "-" + to_string(stop);
    if (sam_hdr_add_line(out_header, "CO", co_value.c_str(), nullptr) != 0) {
        throw runtime_error("Error while writing 'HD' header line to output CRAM file.");
    }
    // END.

    // BEGIN: Copy SQ header lines from input CRAM files to output header.
    samFile *in_cram_fp = open_cram_with_retry(get<0>(begin(crams)->second).c_str());
    if (!in_cram_fp) {
        throw runtime_error("Error while opening input CRAM file.");
    }
    if (hts_set_fai_filename(in_cram_fp, reference_file.c_str()) != 0) {
        throw runtime_error("Error while opening reference FASTA file.");
    }
    sam_hdr_t *in_header = sam_hdr_read(in_cram_fp);
    if (!in_header) {
        throw runtime_error("Error while reading header from input CRAM file.");
    }
    auto in_header_length = sam_hdr_length(in_header);
    istringstream in_header_stream(sam_hdr_str(in_header));
    while (getline(in_header_stream, line)) {
        if (line.find("@SQ", 0) == 0) {
            if (sam_hdr_add_lines(out_header, line.c_str(), line.size()) != 0) {
                throw runtime_error("Error while cpying 'SQ' lines from input CRAM file header to output CRAM file header.");
            }
        }
        if (in_header_stream.tellg() > in_header_length) { // sanity check
            throw runtime_error("Error while parsing header from input CRAM file.");
        }
    }
    sam_hdr_destroy(in_header);
    sam_close(in_cram_fp);
    // END.

    // BEGIN: Write output header to output CRAM file.
    if (sam_hdr_write(out_cram_fp, out_header) < 0) {
        throw runtime_error("Error while writing header to output CRAM file.");
    }
    // END.

    // BEGIN: Process sample by sample.
    unsigned int counter = 0;
    for (auto &&sample: samples) {
        cout << "Processing sample " << sample.first << " (" << ++counter << " out of " << samples.size() << ")" << endl << flush;
        process_sample(get<0>(crams.at(sample.first)), get<1>(crams.at(sample.first)), reference_file, sample.second, out_cram_fp, out_header, variants, window);
    }
    // END.

    // BEGIN: Close output CRAM file.
    sam_hdr_destroy(out_header);
    sam_close(out_cram_fp);
    // END.

    return 0;
}