#include <iostream>
#include <exception>
#include <string>
#include "aux.h"
#include "synced_bcf_reader.h"
#include "TypeSwitcher.h"
#include <vector>
#include <map>
#include <set>
#include <numeric>
#include <boost/program_options.hpp>
#include <vcf.h>
#include <random>
#include <algorithm>

namespace po = boost::program_options;

using namespace std;
using namespace aux;

struct Variant {
    unsigned int sample_idx;
    unsigned int pos;
    bool is_hom;
    string chrom;
    string ref;
    string alt;
    Variant(const string& chrom, unsigned int pos, const string& ref, const string& alt, unsigned int sample_idx, bool is_hom) : chrom(chrom), pos(pos), ref(ref), alt(alt), sample_idx(sample_idx), is_hom(is_hom) {}
};

int main(int argc, char* argv[]) {
    string input_file("");
    string output_file("");
    string samples_file("");
    string samples("");
    string region("");
    int n_random = 0;
    int random_seed = 0;

    random_device rd;

    po::options_description desc("Options");

    desc.add_options()
            ("help,h", "Produce help message")
            ("in,i", po::value<string>(&input_file)->required(), "Input VCF/BCF file. Must be indexed using tabix.")
            ("samples,s", po::value<string>(&samples_file), "Input file with samples. One sample per line.")
            ("region,r", po::value<string>(&region), "Region to be processed. Must follow <CHR>:<START_BP>-<END_BP> format.")
            ("k,k", po::value<int>(&n_random)->default_value(1), "How many individuals to select at random.")
            ("seed,e", po::value<int>(&random_seed)->default_value(rd()), "Random seed.")
            ("out,o", po::value<string>(&output_file)->default_value("variants.tsv.gz"), "Output file name (compressed using bgzip).")
            ;

    po::variables_map vm;

    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        if (vm.count("help")) {
            cout << "This program for every variant extracts K random heterozygous individuals and K random homozygous individuals. Multi-allelic variants are split into bi-allelic variants."  << endl << endl;
            cout << desc << endl;
            return 0;
        }
        po::notify(vm);

        if (n_random <= 0) {
            throw invalid_argument("the option '-k/--k' must be greater than 0");
        }
    } catch (po::error &e) {
        cout << "Error in command line:" << endl;
        cout << e.what() << endl;
        return 1;
    } catch (invalid_argument &e) {
        cout << "Error in command line:" << endl;
        cout << e.what() << endl;
        return 1;
    }

    try {
        if (!samples_file.empty()) {
            samples = read_samples(samples_file.c_str());
        }

        bcf_srs_t *sr = bcf_sr_init();

        if ((input_file.compare("-") != 0) && (!region.empty())) {
            bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);
        }

        if (!region.empty()) {
            if (bcf_sr_set_regions(sr, region.c_str(), 0) < 0) {
                throw runtime_error("Error while subsetting region!");
            }
        }

        if (bcf_sr_add_reader(sr, input_file.c_str()) <= 0) {
            throw runtime_error("Error while initializing VCF/BCF reader!");
        }

        bcf_hdr_t* header = bcf_sr_get_header(sr, 0);

        if (!samples.empty()) {
            if (bcf_hdr_set_samples(header, samples.c_str(), 0) != 0) {
                throw runtime_error("Error while subsetting samples! Check if all specified samples exist in BCF/VCF.");
            }
        }

        int gt_id = bcf_hdr_id2int(header, BCF_DT_ID, "GT");

        if (bcf_hdr_idinfo_exists(header, BCF_HL_FMT, gt_id) == 0) {
            throw runtime_error("GT field was not found!");
        }

        int pl_id = bcf_hdr_id2int(header, BCF_DT_ID, "PL");
        bool pl_exists = (bcf_hdr_idinfo_exists(header, BCF_HL_FMT, pl_id) != 0);

        // don't need pl_exists or move it to options (many times PL is defined in header but not present in variant lines and this cause htslib to crash)
        pl_exists = false;

        int gt_index;
        int pl_index;
        TypeSwitcher gt_switcher;
        TypeSwitcher pl_switcher;
        vector<int32_t> gt_values;
        vector<int32_t> pl_values;

        int allele = 0;
        unsigned int n_chromosomes;
        unsigned int ac_sample, ac_total;
        set<unsigned int> alleles_sample;
        vector<unsigned int> ac;
        map<unsigned int, vector<int>> hom_samples;
        map<unsigned int, vector<int>> het_samples;

        map<string, vector<Variant>> selected_samples;

        mt19937 number_generator(random_seed);

        while (bcf_sr_next_line(sr) > 0) {
            bcf1_t* rec = bcf_sr_get_line(sr, 0);

            if ((sr->streaming == 0) && (rec->pos < sr->regions->start)) {
                continue;
            }

            if ((rec->unpacked & BCF_UN_FMT) == 0) {
                bcf_unpack(rec, BCF_UN_FMT);
            }

            if ((rec->unpacked & BCF_UN_FLT) == 0) {
                bcf_unpack(rec, BCF_UN_FLT);
            }

            gt_index = -1;
            pl_index = -1;

            if (pl_exists) {
                for (int i = 0; i < rec->n_fmt; ++i) {
                    if (rec->d.fmt[i].id == gt_id) {
                        gt_index = i;
                    } else if (rec->d.fmt[i].id == pl_id) {
                        pl_index = i;
                    }
                }

                if (pl_index == -1) {
                    throw runtime_error("PL field was not found in FORMAT!");
                }

                pl_switcher.init(&rec->d.fmt[pl_index]);
            } else {
                for (int i = 0; i < rec->n_fmt; ++i) {
                    if (rec->d.fmt[i].id == gt_id) {
                        gt_index = i;
                        break;
                    }
                }
            }

            if (gt_index == -1) {
                throw runtime_error("GT field was not found in FORMAT!");
            }

            gt_switcher.init(&rec->d.fmt[gt_index]);

            fill(ac.begin(), ac.end(), 0u);
            if (rec->n_allele > ac.size()) {
                ac.resize(rec->n_allele, 0u);
            }
            hom_samples.clear();
            het_samples.clear();

            for (int i = 0; i < rec->n_sample; ++i) {
                if (pl_index >= 0) {
                    (pl_switcher.*(pl_switcher.read))(pl_values);
                    if (accumulate(pl_values.begin(), pl_values.end(), 0) == 0) { // likelihoods for all genotypes are zero
                        continue;
                    }
                }
                (gt_switcher.*(gt_switcher.read))(gt_values);

                n_chromosomes = 0u;
                ac_sample = 0u;
                alleles_sample.clear();
                for (auto&& v : gt_values) {
                    ++n_chromosomes;
                    if (!bcf_gt_is_missing(v)) {
                        allele = bcf_gt_allele(v);
                        ++ac_sample;
                        ac[allele] += 1u;
                        alleles_sample.insert(allele);
                    }
                }

                if ((ac_sample == 0u) || (ac_sample != n_chromosomes)) { // all or some alleles are missing for this sample,
                    continue;
                }

                if (alleles_sample.size() == 1) { // homozygous for some allele
                    if (alleles_sample.erase(0) == 0) { // homozygous for ALT
                        allele = *alleles_sample.begin();
                        auto it = hom_samples.find(allele);
                        if (it != hom_samples.end()) {
                            it->second.push_back(i);
                        } else {
                            hom_samples.emplace(allele, vector<int>({i}));
                        }
                    }
                } else if (alleles_sample.size() > 1) {
                    if (alleles_sample.erase(0) != 0) { // heterozygous with 1 REF and 1 ALT
                        allele =  *alleles_sample.begin();
                        auto it = het_samples.find(allele);
                        if (it != het_samples.end()) {
                            it->second.push_back(i);
                        } else {
                            het_samples.emplace(allele, vector<int>({i}));
                        }
                    }
                }
            }

            ac_total = ac[1];
            for (int i = 2; i < rec->n_allele; ++i) {
                ac_total += ac[i];
            }
            if (ac_total == 0) {
                continue;
            }

            for (int i = 1; i < rec->n_allele; ++i) {
                auto hom_it = hom_samples.find(i);
                if (hom_it != hom_samples.end()) {
                    shuffle(hom_it->second.begin(), hom_it->second.end(), number_generator);
                    for (unsigned int j = 0u; (j < hom_it->second.size() && (j < n_random)); ++j) {
                        selected_samples.emplace(header->samples[hom_it->second[j]], vector<Variant>()).first->second.emplace_back(bcf_seqname(header, rec), rec->pos + 1, rec->d.allele[0], rec->d.allele[i], j + 1, true);
                    }
                }
                auto het_it = het_samples.find(i);
                if (het_it != het_samples.end()) {
                    shuffle(het_it->second.begin(), het_it->second.end(), number_generator);
                    for (unsigned int j = 0u; (j < het_it->second.size()) && (j < n_random); ++j) {
                        selected_samples.emplace(header->samples[het_it->second[j]], vector<Variant>()).first->second.emplace_back(bcf_seqname(header, rec), rec->pos + 1, rec->d.allele[0], rec->d.allele[i], j + 1, false);
                    }
                }
            }
        }

        BGZF *ofp = bgzf_open(output_file.c_str(), "w");
        if (!ofp) {
            throw runtime_error("Error while opening output file!");
        }
        write(ofp, "#RANDOM_SEED=%d\n", random_seed);
        write(ofp, "#MAX_RANDOM_HOM_HETS=%d\n", n_random);
        if (!samples.empty()) {
            write(ofp, "#SAMPLES_USED=%d\n", count(samples.begin(), samples.end(), ',') + 1);
        } else {
            write(ofp, "#SAMPLES_USED=ALL\n");
        }
        write(ofp, "#SAMPLE\tCHROM\tPOS\tREF\tALT\tIS_HOM\tSAMPLE_INDEX\n");
        for (auto &&selected_sample: selected_samples) {
            for (auto &&variant : selected_sample.second) {
                write(ofp, "%s\t%s\t%lu\t%s\t%s\t%d\t%d\n", selected_sample.first.c_str(), variant.chrom.c_str(), variant.pos, variant.ref.c_str(), variant.alt.c_str(), variant.is_hom, variant.sample_idx);
            }
        }
        if (bgzf_close(ofp) != 0) {
            throw runtime_error("Error while closing output file!");
        }
        bcf_sr_destroy(sr);
    } catch (exception &e) {
        cout << "Error: " << endl;
        cout << e.what() << endl;
        return 1;
    }

    return 0;
}
