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
#include "bgzf.h"
#include "htslib/sam.h"

namespace po = boost::program_options;

using namespace std;

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

vector<string> variants_file_header {
        "#SAMPLE", "CHROM", "POS", "REF", "ALT", "IS_HOM", "SAMPLE_INDEX"
};

struct Region {
    bool is_hom;
    unsigned int sample_idx;
    unsigned int variant_pos;
    string variant_ref;
    string variant_alt;
    string chrom;
    Region(bool is_hom, unsigned int sample_idx, unsigned int variant_pos, const string& variant_ref, const string& variant_alt, const string& chrom): is_hom(is_hom), sample_idx(sample_idx), variant_pos(variant_pos), variant_ref(variant_ref), variant_alt(variant_alt), chrom(chrom) {}
    bool operator<(const Region& region) const {
        int chrom_cmp = chrom.compare(region.chrom);
        return chrom_cmp == 0 ? variant_pos < region.variant_pos : chrom_cmp < 0;
    }
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

int main(int argc, char* argv[]) {
    string cram_file("");
    string crai_file("");
    string variants_file("");
    string reference_file("");
    string output_cram_file("");
    unsigned int window = 0u;

    po::options_description desc("Options");

    desc.add_options()
            ("help,h", "Produce help message")
            ("cram,c", po::value<string>(&cram_file)->required(), "Single sample CRAM file.")
            ("crai,i", po::value<string>(&crai_file)->required(), "CRAM index file (CRAI).")
            ("variants,v", po::value<string>(&variants_file)->required(), "Compressed (gzip) file with a list of variants to extract from sample CRAM.")
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

    // BEGIN: Read variants for sample.
    vector<Region> regions;

    BGZF *if_variants = bgzf_open(variants_file.c_str(), "r");
    if (!if_variants) {
        throw runtime_error("Error while opening variants file!");
    }
    kstring_t kline = KS_INITIALIZE;
    bool found_header = false;
    vector<string> tokens;
    auto separator = regex("[ \t]");
    while (bgzf_getline(if_variants, '\n', &kline) >= 0) {
        tokens.clear();
        copy(cregex_token_iterator(kline.s, kline.s + kline.l, separator, -1), cregex_token_iterator(), back_inserter(tokens));
        if ((tokens.size() == 0) || (tokens[0].empty()) || (tokens[0][0] == '#')) { // header lines
            if ((tokens.size() > 0) && (tokens[0].compare(variants_file_header[0]) == 0)) {
                found_header = true;
                if (tokens.size() != variants_file_header.size()) {
                    throw runtime_error("Incorrect header in variants file!");
                }
                for (unsigned int i = 0u; i < tokens.size(); ++i) {
                    if (tokens[i].compare(variants_file_header[i]) != 0) {
                        throw runtime_error("Incorrect header in variants file!");
                    }
                }
            }
            continue;
        }
        if (!found_header) {
            throw runtime_error("No header in variants file!");
        }
        if (tokens.size() != variants_file_header.size()) {
            throw runtime_error("Incorrect number of columns in variants file.");
        }
        // 0 - "#SAMPLE", 1 - "CHROM", 2 - "POS", 3 -  "REF", 4 - "ALT", 5 - "IS_HOM", 6 - "SAMPLE_INDEX"
        // cout << (tokens[5].compare("1") == 0) << " " << stoul(tokens[6]) << " " << stoul(tokens[2]) << " " << tokens[3] << " " << tokens[4] << " " << tokens[1] << endl;

        regions.emplace_back(tokens[5].compare("1") == 0, stoul(tokens[6]), stoul(tokens[2]), tokens[3], tokens[4], tokens[1]);
    }
    free(ks_release(&kline));
    if (bgzf_close(if_variants) != 0) {
        throw runtime_error("Error while closing variants file!");
    }
    // END.

    // BEGIN: Open input CRAM.
    samFile *in_cram_fp = open_cram_with_retry(cram_file.c_str());
    if (!in_cram_fp) {
        throw runtime_error("Error while opening input CRAM file.");
    }
    if (hts_set_fai_filename(in_cram_fp, reference_file.c_str()) != 0) {
        throw runtime_error("Error while opening reference FASTA file.");
    }
    bam_hdr_t *in_header = sam_hdr_read(in_cram_fp);
    if (!in_header) {
        throw runtime_error("Error while reading header from input CRAM file.");
    }
    hts_idx_t *in_idx = open_index_with_retry(in_cram_fp, cram_file.c_str(), crai_file.c_str());
    if (!in_idx) {
        throw runtime_error("Error while reading input CRAI file.");
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
    // END.

    // BEGIN: Copy SQ header lines from input CRAM files to output header.
    auto in_header_length = sam_hdr_length(in_header);
    istringstream in_header_stream(sam_hdr_str(in_header));
    string line("");
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
    // END.

    // BEGIN: Write output header to output CRAM file.
    if (sam_hdr_write(out_cram_fp, out_header) < 0) {
        throw runtime_error("Error while writing header to output CRAM file.");
    }
    // END.

    // BEGIN: Copy sequences from input CRAM to output CRAM.
    sort(regions.begin(), regions.end());
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
    for (auto &&region: regions) {
        int tid = bam_name2id(in_header, region.chrom.c_str());
        if (tid < 0) {
            throw runtime_error("Error while getting reference id from input CRAM file.");
        }
        start = window > region.variant_pos ? 0u : region.variant_pos - window;
        stop = region.variant_pos + window;
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
            new_qname = to_string(region.variant_pos) + ":" + region.variant_ref + ":" + region.variant_alt + ":" + (region.is_hom ? "0" : "") + to_string(region.sample_idx) + ":" + to_string(qname_idx);
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
    // END.

    // BEGIN: Close open handlers.
    bam_destroy1(b);
    hts_idx_destroy(in_idx);
    sam_hdr_destroy(in_header);
    sam_close(in_cram_fp);
    sam_hdr_destroy(out_header);
    sam_close(out_cram_fp);
    // END.

    return 0;
}