#include <iostream>
#include <fstream>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <algorithm>
#include <vector>
#include <set>
#include "hts.h"
#include "vcf.h"
#include "synced_bcf_reader.h"
#include "Histogram.h"
#include "TypeSwitcher.h"
#include "aux.h"
#include <chrono>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

using namespace std;
using namespace aux;

int main(int argc, char* argv[]) {
	string input_gt_file("");
    string input_dp_file("");
    string samples_file("");
	string output_file("");
	string samples("");

	vector<double> hist_borders = { 0, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51, 56, 61, 66, 71, 76, 81, 86, 91, 96, numeric_limits<double>::max() };

	po::options_description desc("Options");

	desc.add_options()
		("help,h", "Produce help message")
		("genotypes,g", po::value<string>(&input_gt_file)->required(), "Input VCF/BCF file with GT.")
		("depth,d", po::value<string>(&input_dp_file)->required(), "Input VCF/BCF file with DP.")
		("samples,s", po::value<string>(&samples_file), "Input file with samples. One sample per line.")
		("out,o", po::value<string>(&output_file)->required(), "Output file. Compressed using gzip.")
	;

	po::variables_map vm;

	try {
		po::store(po::parse_command_line(argc, argv, desc), vm);
		if (vm.count("help")) {
			cout << "For each variant this program computes: NS, AN, AC, AF, Hom, Het, AVGDP, AVGGQ, histograms of depth (DP) and genotype qualities (GQ)." << endl << endl;
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
            cout << "Reading samples file... " << flush;
            samples = read_samples(samples_file.c_str());
            cout << "Done." << endl;
            cout << "Found " << count(samples.begin(), samples.end(), ',')  + 1 << " sample(s)." << endl; 
        }

        BGZF *ofp = bgzf_open(output_file.c_str(), "w");
        if (!ofp) {
            throw runtime_error("Error while opening output file!");
        }

		bcf_srs_t *gt_sr = bcf_sr_init();
        bcf_srs_t *dp_sr = bcf_sr_init();

		if ((bcf_sr_add_reader(gt_sr, input_gt_file.c_str()) <= 0) || (bcf_sr_add_reader(dp_sr, input_dp_file.c_str()) <= 0)) {
			throw runtime_error("Error while initializing VCF/BCF reader!");
		}

		bcf_hdr_t* gt_sr_header = bcf_sr_get_header(gt_sr, 0);
        bcf_hdr_t* dp_sr_header = bcf_sr_get_header(dp_sr, 0);

		if (!samples.empty()) {
		    cout << "Setting sample(s)... " << flush;
			if ((bcf_hdr_set_samples(gt_sr_header, samples.c_str(), 0) != 0) || (bcf_hdr_set_samples(dp_sr_header, samples.c_str(), 0) != 0)) {
				throw runtime_error("Error while subsetting samples!");
			}
            cout << "Done." << endl;
		}

		int gt_id = bcf_hdr_id2int(gt_sr_header, BCF_DT_ID, "GT");
		int dp_id = bcf_hdr_id2int(dp_sr_header, BCF_DT_ID, "DP");

		if (bcf_hdr_idinfo_exists(gt_sr_header, BCF_HL_FMT, gt_id) == 0) {
			throw runtime_error("GT field was not found!");
		}

		if (bcf_hdr_idinfo_exists(dp_sr_header, BCF_HL_FMT, dp_id) == 0) {
			throw runtime_error("DP field was not found!");
		}

		if (bcf_hdr_nsamples(gt_sr_header) != bcf_hdr_nsamples(dp_sr_header)) {
		    throw runtime_error("Different number of samples in VCF/BCF files.");
		}

		cout << "Using " << bcf_hdr_nsamples(dp_sr_header) << " sample(s)" << endl;

		for (unsigned int i = 0u; i < bcf_hdr_nsamples(gt_sr_header); ++i) {
		    if (strcmp(gt_sr_header->samples[i], dp_sr_header->samples[i]) != 0) {
		        throw runtime_error("Samples or their order don't match in VCF/BCF files.");
		    }
		}

        write(ofp, "##fileformat=VCFv4.2\n");
        for (int i = 0; i < gt_sr_header->nhrec; ++i) {
            if (strcmp(gt_sr_header->hrec[i]->key, "FILTER") == 0) {
                write(ofp, "##%s=<%s=%s", gt_sr_header->hrec[i]->key, gt_sr_header->hrec[i]->keys[0], gt_sr_header->hrec[i]->vals[0]);
                for (int j = 1; j < gt_sr_header->hrec[i]->nkeys; ++j) {
                    if (strcmp(gt_sr_header->hrec[i]->keys[j], "IDX") == 0) {
                        continue;
                    }
                    write(ofp, ",%s=%s", gt_sr_header->hrec[i]->keys[j], gt_sr_header->hrec[i]->vals[j]);
                }
                write(ofp, ">\n");
            }
        }
        write(ofp, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Coverage\">\n");
        write(ofp, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Number of Alleles in Samples with Coverage\">\n");
        write(ofp, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate Allele Counts in Samples with Coverage\">\n");
        write(ofp, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternate Allele Frequencies\">\n");
        write(ofp, "##INFO=<ID=Het,Number=A,Type=Integer,Description=\"Heterozygous Counts\">\n");
        write(ofp, "##INFO=<ID=Hom,Number=A,Type=Integer,Description=\"Homozygous Alternate Counts\">\n");
        write(ofp, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth at site\">\n");
		  write(ofp, "##INFO=<ID=DP_R,Number=R,Type=Integer,Description=\"Total depth at site accross samples carrying allele\">\n");
        write(ofp, "##INFO=<ID=AVGDP,Number=1,Type=Float,Description=\"Average depth per sample\">\n");
        write(ofp, "##INFO=<ID=AVGDP_R,Number=R,Type=Float,Description=\"Average depth per sample carrying allele\">\n");
        write(ofp, "##INFO=<ID=DP_HIST,Number=1,Type=String,Description=\"Histogram of DP across all samples; Mids: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5\">\n");
		write(ofp, "##INFO=<ID=DP_HIST_R,Number=R,Type=String,Description=\"Histograms of DP across samples carrying allele; Mids: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5\">\n");
		write(ofp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

		int gt_index, dp_index;
		TypeSwitcher gt_switcher, dp_switcher;
		unsigned int ns, an, ac_sample, ac_total;
        set<unsigned int> hom_sample;
		vector<unsigned int> ac, hom, het;
        int allele = 0;
		vector<int32_t> gt_values, dp_values;
		vector<Histogram> dp_histograms;
		set<int> unique_alleles;

        cout << "Processing... " << flush;
		while ((bcf_sr_next_line(gt_sr) > 0) && (bcf_sr_next_line(dp_sr) > 0)) {
         //auto start = std::chrono::system_clock::now();
         bcf1_t* gt_rec = bcf_sr_get_line(gt_sr, 0);
            bcf1_t* dp_rec = bcf_sr_get_line(dp_sr, 0);

            if (gt_rec->pos < dp_rec->pos) {
                while (bcf_sr_next_line(gt_sr) > 0) {
                    gt_rec = bcf_sr_get_line(gt_sr, 0);
                    if (gt_rec->pos == dp_rec->pos) {
                        if (gt_rec->n_allele == dp_rec->n_allele) {
                            bool alleles_match = true;
                            for (int i = 0; i < gt_rec->n_allele; ++i) {
                                if (strcmp(gt_rec->d.allele[i], dp_rec->d.allele[i]) != 0) {
                                    alleles_match = false;
                                    break;
                                }
                            }
                            if (alleles_match) {
                                break;
                            }
                        }
                    } else if (gt_rec->pos > dp_rec->pos) {
                        throw runtime_error("Incompatible VCF/BCF files.");
                    }
                }
            }

            if (gt_rec->pos != dp_rec->pos) {
                throw runtime_error("Found incompatible records: positions don't match.");
            }
            if (gt_rec->n_allele != dp_rec->n_allele) {
                throw runtime_error("Found incompatible records: number of alleles don't match.");
            }
            for (int i = 0; i < gt_rec->n_allele; ++i) {
                if (strcmp(gt_rec->d.allele[i], dp_rec->d.allele[i]) != 0) {
                    throw runtime_error("Found incompatible records: alleles don't match.");
                }
            }

            if ((gt_rec->unpacked & BCF_UN_FMT) == 0) {
                bcf_unpack(gt_rec, BCF_UN_FMT);
            }

            if ((gt_rec->unpacked & BCF_UN_FLT) == 0) {
                bcf_unpack(gt_rec, BCF_UN_FLT);
            }

            if ((dp_rec->unpacked & BCF_UN_FMT) == 0) {
                bcf_unpack(dp_rec, BCF_UN_FMT);
            }

            if ((dp_rec->unpacked & BCF_UN_FLT) == 0) {
                bcf_unpack(dp_rec, BCF_UN_FLT);
            }

            gt_index = -1;
            for (int i = 0; i < gt_rec->n_fmt; ++i) {
                if (gt_rec->d.fmt[i].id == gt_id) {
                    gt_index = i;
                    break;
                }
            }
            if (gt_index == -1) {
                throw runtime_error("GT field was not found in FORMAT!");
            }
            gt_switcher.init(&gt_rec->d.fmt[gt_index]);

            dp_index = -1;
            for (int i = 0; i < dp_rec->n_fmt; ++i) {
                if (dp_rec->d.fmt[i].id == dp_id) {
                    dp_index = i;
                    break;
                }
            }
			if (dp_index == -1) {
                cerr << "[warning] DP field missing for " << bcf_seqname(dp_sr_header, dp_rec) << ":" << dp_rec->pos + 1 << ":" << dp_rec->d.allele[0];
                for (int i = 1; i < dp_rec->n_allele; ++i) {
                    cerr << "/" << dp_rec->d.allele[i];
                }
                cerr << endl;
			} else {
                dp_switcher.init(&dp_rec->d.fmt[dp_index]);
			}

            ns = 0u;
            an = 0u;
			fill(ac.begin(), ac.end(), 0u); // cleanup allele, hom and het counts
            fill(hom.begin(), hom.end(), 0u);
            fill(het.begin(), het.end(), 0u);
			if (gt_rec->n_allele > ac.size()) { // append additional allele, hom and het counts if needed
				ac.resize(gt_rec->n_allele, 0u);
                hom.resize(gt_rec->n_allele, 0u);
                het.resize(gt_rec->n_allele, 0u);
			}

			for (auto&& h : dp_histograms) { // cleanup DP histograms and reuse them
				h.clear();
			}
			for (int i = dp_histograms.size(); i < dp_rec->n_allele; ++i) { // append additional DP histograms if needed
				dp_histograms.emplace_back(Histogram(hist_borders));
			}

            Histogram dp_histogram(hist_borders);

			for (int i = 0; i < gt_rec->n_sample; ++i) {
				(gt_switcher.*(gt_switcher.read))(gt_values);
				if (dp_index != -1) (dp_switcher.*(dp_switcher.read))(dp_values);

				ac_sample = 0u;
                hom_sample.clear();
				unique_alleles.clear();
				for (auto&& v : gt_values) {
					if (!bcf_gt_is_missing(v)) {
					    allele = bcf_gt_allele(v);
						unique_alleles.insert(allele);
						ac[allele] += 1u;
                        ++ac_sample;
                        hom_sample.insert(allele);
					}
				}

                if (ac_sample == 0u) { // all alleles missing for this sample
                    continue;
                }

                ++ns;
                an += ac_sample;

                if (hom_sample.size() == 1) { // homozygous for some allele
                    hom[*hom_sample.begin()] += 1;
                } else if (hom_sample.size() > 1) {
                    if (hom_sample.erase(0) != 0) { // heterozygous with 1 REF and 1 ALT
                        het[*hom_sample.begin()] += 1;
                    }
                }

                if (dp_index != -1 && dp_values[0] != bcf_int32_missing) {
                    dp_histogram.add((int)dp_values[0]);
                }

                for (auto&& allele: unique_alleles) {
				    if (dp_index != -1 && dp_values[0] != bcf_int32_missing) {
					    dp_histograms[allele].add((int)dp_values[0]);
				    }
				}
			}

			ac_total = ac[1];
			for (int i = 2; i < gt_rec->n_allele; ++i) {
				ac_total += ac[i];
			}
			//if (ac_total == 0) {
			//	continue;
			//}

         //auto end = std::chrono::system_clock::now();
         //std::chrono::duration<double> elapsed = end - start;
         //std::cout << "Segment load elapsed time: " << elapsed.count() << " s\n";

			write(ofp, "%s\t%lu\t%s\t%s\t%s", bcf_seqname(gt_sr_header, gt_rec), gt_rec->pos + 1, gt_rec->d.id, gt_rec->d.allele[0], gt_rec->d.allele[1]);
			for (int i = 2; i < gt_rec->n_allele; ++i) {
				write(ofp, ",%s", gt_rec->d.allele[i]);
			}
            if (isnan(gt_rec->qual)) {
                write(ofp, "\t.");
            } else {
                write(ofp, "\t%g", gt_rec->qual);
            }
            if (gt_rec->d.n_flt <= 0) {
                write(ofp, "\t.");
            } else {
                write(ofp, "\t%s", bcf_hdr_int2id(gt_sr_header, BCF_DT_ID, gt_rec->d.flt[0]));
                for (int i = 1; i < gt_rec->d.n_flt; ++i) {
                    write(ofp, ";%s", bcf_hdr_int2id(gt_sr_header, BCF_DT_ID, gt_rec->d.flt[i]));
                }
            }
            write(ofp, "\tNS=%d;AN=%d", ns, an);
            write(ofp, ";AC=%d", ac[1]);
            for (int i = 2; i < gt_rec->n_allele; ++i) {
                write(ofp, ",%d", ac[i]);
            }
            write(ofp, ";AF=%g", an > 0 ? (ac[1] / (double)an) : 0.0);
            for (int i = 2; i < gt_rec->n_allele; ++i) {
                write(ofp, ",%g", an > 0 ? (ac[i] / (double)an) : 0.0);
            }
            write(ofp, ";Het=%d", het[1]);
            for (int i = 2; i < gt_rec->n_allele; ++i) {
                write(ofp, ",%d", het[i]);
            }
            write(ofp, ";Hom=%d", hom[1]);
            for (int i = 2; i < gt_rec->n_allele; ++i) {
                write(ofp, ",%d", hom[i]);
            }
            if (dp_index != -1) {
               write(ofp, ";DP=%d", (long long int)dp_histogram.get_total());
               write(ofp, ";DP_R=%d", (long long int)dp_histograms[0].get_total());
               for (int i = 1; i < gt_rec->n_allele; ++i) {
                  if (dp_histograms[i].get_n() > 0) {
                     write(ofp, ",%d", (long long int)dp_histograms[i].get_total());
                  } else {
                     write(ofp, ",.");
                  }
               }
               write(ofp, ";AVGDP=%g", dp_histogram.get_average());
               write(ofp, ";AVGDP_R=%g", dp_histograms[0].get_average());
               for (int i = 1; i < gt_rec->n_allele; ++i) {
                  if (dp_histograms[i].get_n() > 0) {
                     write(ofp, ",%g", dp_histograms[i].get_average());
                  } else {
                     write(ofp, ",.");
                  }
               }
            }
            if (dp_index != -1) {
               write(ofp, ";DP_HIST=%s", dp_histogram.get_text());
               write(ofp, ";DP_HIST_R=%s", dp_histograms[0].get_text());
               for (int i = 1; i < gt_rec->n_allele; ++i) {
                  write(ofp, ",%s", dp_histograms[i].get_text());
               }
            }
			write(ofp, "\n");
		}
        cout << "Done." << endl;

		bcf_sr_destroy(gt_sr);
        bcf_sr_destroy(dp_sr);

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
