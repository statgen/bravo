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
#include "Histogram.h"
#include "TypeSwitcher.h"
#include "aux.h"

#include <boost/program_options.hpp>

namespace po = boost::program_options;

using namespace std;
using namespace aux;

int main(int argc, char* argv[]) {
	string input_file("");
	string output_file("");
	string samples_file("");
	string samples("");
	string region("");

	vector<double> hist_borders = { 0, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51, 56, 61, 66, 71, 76, 81, 86, 91, 96, numeric_limits<double>::max() };

	po::options_description desc("Options");

	desc.add_options()
		("help,h", "Produce help message")
		("in,i", po::value<string>(&input_file)->required(), "Input VCF/BCF file. Must be indexed using tabix.")
		("samples,s", po::value<string>(&samples_file), "Input file with samples. One sample per line.")
		("region,r", po::value<string>(&region), "Region to be processed. Must follow <CHR>:<START_BP>-<END_BP> format.")
		("out,o", po::value<string>(&output_file)->required(), "Output file. Compressed using gzip.")
	;

	po::variables_map vm;

	try {
		po::store(po::parse_command_line(argc, argv, desc), vm);
		if (vm.count("help")) {
			cout << "This program computes histograms of depth (DP) and genotype qualities (GQ) for each variant." << endl << endl;
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
        if (!samples_file.empty()) {
            samples = read_samples(samples_file.c_str());
        }

        BGZF *ofp = bgzf_open(output_file.c_str(), "w");
        if (!ofp) {
            throw runtime_error("Error while opening output file!");
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
				throw runtime_error("Error while subsetting samples!");
			}
		}

		int gt_id = bcf_hdr_id2int(header, BCF_DT_ID, "GT");
		int dp_id = bcf_hdr_id2int(header, BCF_DT_ID, "DP");
		int gq_id = bcf_hdr_id2int(header, BCF_DT_ID, "GQ");

		if (bcf_hdr_idinfo_exists(header, BCF_HL_FMT, gt_id) == 0) {
			throw runtime_error("GT field was not found!");
		}

		if (bcf_hdr_idinfo_exists(header, BCF_HL_FMT, dp_id) == 0) {
			throw runtime_error("DP field was not found!");
		}

		if (bcf_hdr_idinfo_exists(header, BCF_HL_FMT, gq_id) == 0) {
			throw runtime_error("GQ field was not found!");
		}

		write(ofp, "##fileformat=VCFv4.2\n");

        write(ofp, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth at site\">\n");
		write(ofp, "##INFO=<ID=AVGDP,Number=1,Type=Float,Description=\"Average depth per sample\">\n");
        write(ofp, "##INFO=<ID=AVGDP_R,Number=R,Type=Float,Description=\"Average depth per sample carrying allele\">\n");
        write(ofp, "##INFO=<ID=AVGGQ,Number=1,Type=Float,Description=\"Average genotype quality per sample\">\n");
        write(ofp, "##INFO=<ID=AVGGQ_R,Number=R,Type=Float,Description=\"Average genotype quality per sample carrying allele\">\n");
        write(ofp, "##INFO=<ID=DP_HIST,Number=1,Type=String,Description=\"Histogram of DP across all samples; Mids: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5\">\n");
		write(ofp, "##INFO=<ID=DP_HIST_R,Number=R,Type=String,Description=\"Histograms of DP across samples carrying allele; Mids: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5\">\n");
		write(ofp, "##INFO=<ID=GQ_HIST,Number=1,Type=String,Description=\"Histogram of GQ across all samples; Mids: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5\">\n");
        write(ofp, "##INFO=<ID=GQ_HIST_R,Number=R,Type=String,Description=\"Histograms of GQ across samples carrying allele; Mids: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5\">\n");
		write(ofp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

		int gt_index, dp_index, gq_index;
		TypeSwitcher gt_switcher, dp_switcher, gq_switcher;
		unsigned int ac_sample, ac_total;
		vector<unsigned int> ac;
        int allele = 0;
		vector<int32_t> gt_values, dp_values, gq_values;
		vector<Histogram> dp_histograms, gq_histograms;
		set<int> unique_alleles;

		while (bcf_sr_next_line(sr) > 0) {
			bcf1_t* rec = bcf_sr_get_line(sr, 0);

            if ((sr->streaming == 0) && (rec->pos < sr->regions->start)) {
                continue;
            }

			if (rec->unpacked == 1) {
				bcf_unpack(rec, BCF_UN_FMT);
			}

			gt_index = -1;
			dp_index = -1;
			gq_index = -1;

			for (int i = 0; i < rec->n_fmt; ++i) {
				if (rec->d.fmt[i].id == gt_id) {
					gt_index = i;
				} else if (rec->d.fmt[i].id == dp_id) {
					dp_index = i;
				} else if (rec->d.fmt[i].id == gq_id) {
					gq_index = i;
				}
			}

			if (gt_index == -1) {
				throw runtime_error("GT field was not found in FORMAT!");
			}

			if (dp_index == -1) {
				throw runtime_error("DP field was not found in FORMAT!");
			}

			if (gq_index == -1) {
				throw runtime_error("GQ field was not found in FORMAT!");
			}

			gt_switcher.init(&rec->d.fmt[gt_index]);
			dp_switcher.init(&rec->d.fmt[dp_index]);
			gq_switcher.init(&rec->d.fmt[gq_index]);

			fill(ac.begin(), ac.end(), 0u); // cleanup allele counts
			if (rec->n_allele > ac.size()) { // append additional allele counts if needed
				ac.resize(rec->n_allele, 0u);
			}

			for (auto&& h : dp_histograms) { // cleanup DP histograms and reuse them
				h.clear();
			}
			for (int i = dp_histograms.size(); i < rec->n_allele; ++i) { // append additional DP histograms if needed
				dp_histograms.emplace_back(Histogram(hist_borders));
			}

			for (auto &&h: gq_histograms) { // cleanup DP histohgrams and reuse them
				h.clear();
			}
			for (int i = gq_histograms.size(); i < rec->n_allele; ++i) { // append additional GQ histograms if needed
				gq_histograms.emplace_back(Histogram(hist_borders));
			}

            Histogram dp_histogram(hist_borders), gq_histogram(hist_borders);

			for (int i = 0; i < rec->n_sample; ++i) {
				(gt_switcher.*(gt_switcher.read))(gt_values);
				(dp_switcher.*(dp_switcher.read))(dp_values);
				(gq_switcher.*(gq_switcher.read))(gq_values);

                ac_sample = 0u;
				unique_alleles.clear();
				for (auto&& v : gt_values) {
					if (!bcf_gt_is_missing(v)) {
					    allele = bcf_gt_allele(v);
						unique_alleles.insert(allele);
						ac[allele] += 1u;
                        ++ac_sample;
					}
				}

                if (ac_sample == 0u) { // all alleles missing for this sample
                    continue;
                }

				if (dp_values[0] != bcf_int32_missing) {
					dp_histogram.add((int)dp_values[0]);
				}

				if (gq_values[0] != bcf_int32_missing) {
				    gq_histogram.add((int)gq_values[0]);
				}

				for (auto&& allele: unique_alleles) {
					if (dp_values[0] != bcf_int32_missing) {
						dp_histograms[allele].add((int)dp_values[0]);
					}
					if (gq_values[0] != bcf_int32_missing) {
						gq_histograms[allele].add((int)gq_values[0]);
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

			write(ofp, "%s\t%lu\t%s\t%s\t%s", bcf_seqname(header, rec), rec->pos + 1, rec->d.id, rec->d.allele[0], rec->d.allele[1]);
			for (int i = 2; i < rec->n_allele; ++i) {
				write(ofp, ",%s", rec->d.allele[i]);
			}
			write(ofp, "\t.\t.\t");

            write(ofp, "DP=%d", (long long int)dp_histogram.get_total());
            write(ofp, ";AVGDP=%g", dp_histogram.get_average());
            write(ofp, ";AVGDP_R=%g", dp_histograms[0].get_average());
            for (int i = 1; i < rec->n_allele; ++i) {
                write(ofp, ",%g", dp_histograms[i].get_average());
            }
            write(ofp, ";AVGGQ=%g", gq_histogram.get_average());
            write(ofp, ";AVGGQ_R=%g", gq_histograms[0].get_average());
            for (int i = 1; i < rec->n_allele; ++i) {
                write(ofp, ",%g", gq_histograms[i].get_average());
            }
            write(ofp, ";DP_HIST=%s", dp_histogram.get_text());
            write(ofp, ";GQ_HIST=%s", gq_histogram.get_text());
            write(ofp, ";DP_HIST_R=%s", dp_histograms[0].get_text());
			for (int i = 1; i < rec->n_allele; ++i) {
				write(ofp, ",%s", dp_histograms[i].get_text());
			}
			write(ofp, ";GQ_HIST_R=%s", gq_histograms[0].get_text());
			for (int i = 1; i < rec->n_allele; ++i) {
				write(ofp, ",%s", gq_histograms[i].get_text());
			}
			write(ofp, "\n");
		}

		bcf_sr_destroy(sr);

        if (bgzf_close(ofp) != 0) {
            throw runtime_error("Error while closing output file!");
        }
	} catch (exception &e) {
		cout << "Error: " << endl;
		cout << e.what() << endl;
      return 1;
   }

	return 0;
}
