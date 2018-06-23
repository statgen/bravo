#include <iostream>
#include <fstream>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <vector>
#include <set>
#include <cmath>
#include <numeric>
#include "hts.h"
#include "vcf.h"
#include "synced_bcf_reader.h"
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
	string label("");
	unsigned int max_allele_length = numeric_limits<unsigned int>::max();

	po::options_description desc("Options");

	desc.add_options()
		("help,h", "Produce help message")
		("in,i", po::value<string>(&input_file)->required(), "Input VCF/BCF file. Must be indexed using tabix.")
		("samples,s", po::value<string>(&samples_file), "Input file with samples. One sample per line.")
		("region,r", po::value<string>(&region), "Region to be processed. Must follow <CHR>:<START_BP>-<END_BP> format.")
		("label,l", po::value<string>(&label), "Study/population label (for dbGaP submission VCF).")
		("bp,b", po::value<unsigned int>(&max_allele_length), "Maximal length of deletions/insertions in base pairs (including the preceding anchor base).")
		("out,o", po::value<string>(&output_file)->required(), "Output file. Compressed using gzip.")
	;

	po::variables_map vm;

	try {
		po::store(po::parse_command_line(argc, argv, desc), vm);
		if (vm.count("help")) {
			cout << "This program computes NS, AN, AC, AF, Hom, Het values for each variant. Output is suitable for a dbGaP submission."  << endl << endl;
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
				throw runtime_error("Error while subsetting samples! Check if all specified samples exist in BCF/VCF.");
			}
		}

		int gt_id = bcf_hdr_id2int(header, BCF_DT_ID, "GT");

		if (bcf_hdr_idinfo_exists(header, BCF_HL_FMT, gt_id) == 0) {
			throw runtime_error("GT field was not found!");
		}

		int pl_id = bcf_hdr_id2int(header, BCF_DT_ID, "PL");
		bool pl_exists = (bcf_hdr_idinfo_exists(header, BCF_HL_FMT, pl_id) != 0);

		write(ofp, "##fileformat=VCFv4.2\n");
		for (int i = 0; i < header->nhrec; ++i) {
			if (strcmp(header->hrec[i]->key, "FILTER") == 0) {
				write(ofp, "##%s=<%s=%s", header->hrec[i]->key, header->hrec[i]->keys[0], header->hrec[i]->vals[0]);
				for (int j = 1; j < header->hrec[i]->nkeys; ++j) {
					if (strcmp(header->hrec[i]->keys[j], "IDX") == 0) {
						continue;
					}
					write(ofp, ",%s=%s", header->hrec[i]->keys[j], header->hrec[i]->vals[j]);
				}
				write(ofp, ">\n");
			}
		}
		if (!label.empty()) {
 			write(ofp, "##FORMAT=<ID=NA,Number=1,Type=Integer,Description=\"Number of alleles for the population.\"\n");
 			write(ofp, "##FORMAT=<ID=FRQ,Number=.,Type=Float,Description=\"Frequency of each alternate allele.\"\n");
	   	}
		write(ofp, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Coverage\">\n");
		write(ofp, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Number of Alleles in Samples with Coverage\">\n");
		write(ofp, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate Allele Counts in Samples with Coverage\">\n");
		write(ofp, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternate Allele Frequencies\">\n");
		write(ofp, "##INFO=<ID=Het,Number=A,Type=Integer,Description=\"Heterozygous Counts\">\n");
		write(ofp, "##INFO=<ID=Hom,Number=A,Type=Integer,Description=\"Homozygous Alternate Counts\">\n");
		if (!label.empty()) {
			write(ofp, "##INFO=<ID=VRT,Number=1,Type=Integer,Description=\"Variation type,1 - SNV: single nucleotide variation,2 - DIV: deletion/insertion variation,3 - HETEROZYGOUS: variable, but undefined at nucleotide level,4 - STR: short tandem repeat (microsatellite) variation, 5 - NAMED: insertion/deletion variation of named repetitive element,6 - NO VARIATON: sequence scanned for variation, but none observed,7 - MIXED: cluster contains submissions from 2 or more allelic classes (not used),8 - MNV: multiple nucleotide variation with alleles of common length greater than 1,9 - Exception\">\n");
		}

		if (!label.empty()) {
			write(ofp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" , label.c_str());
		} else {
			write(ofp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
		}

		int gt_index;
		int pl_index;
		TypeSwitcher gt_switcher;
		TypeSwitcher pl_switcher;
		vector<int32_t> gt_values;
		vector<int32_t> pl_values;
		unsigned int ns, an, ac_total;
		unsigned int ac_sample;
		set<unsigned int> hom_sample;
		vector<unsigned int> ac;
		vector<unsigned int> hom;
		vector<unsigned int> het;
		int allele = 0;
		int vrt = 1;
		set<unsigned int> allele_lengths;

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

			allele_lengths.clear();
			for (int i = 0; i < rec->n_allele; ++i) {
				allele_lengths.insert(strlen(rec->d.allele[i]));
			}

			if (*allele_lengths.rbegin() > max_allele_length) {
				continue;
			}

			if (allele_lengths.size() == 1) {
				if (*allele_lengths.begin() == 1) {
					vrt = 1;
				} else {
					vrt = 8;
				}
			} else {
				vrt = 2;
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

			ns = 0u;
			an = 0u;

			fill(ac.begin(), ac.end(), 0u);
			fill(hom.begin(), hom.end(), 0u);
			fill(het.begin(), het.end(), 0u);
			if (rec->n_allele > ac.size()) {
				ac.resize(rec->n_allele, 0u);
				hom.resize(rec->n_allele, 0u);
				het.resize(rec->n_allele, 0u);
			}

			for (int i = 0; i < rec->n_sample; ++i) {
				if (pl_index >= 0) {
					(pl_switcher.*(pl_switcher.read))(pl_values);
					if (accumulate(pl_values.begin(), pl_values.end(), 0) == 0) { // likelihoods for all genotypes are zero
						continue;
					}
				}

				(gt_switcher.*(gt_switcher.read))(gt_values);

				ac_sample = 0u;
				hom_sample.clear();
				for (auto&& v : gt_values) {
					if (!bcf_gt_is_missing(v)) {
						allele = bcf_gt_allele(v);
						++ac_sample;
						ac[allele] += 1u;
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
			if (isnan(rec->qual)) {
				write(ofp, "\t.");
			} else {
				write(ofp, "\t%g", rec->qual);
			}
			if (rec->d.n_flt <= 0) {
				write(ofp, "\t.");
			} else {
				write(ofp, "\t%s", bcf_hdr_int2id(header, BCF_DT_ID, rec->d.flt[0]));
				for (int i = 1; i < rec->d.n_flt; ++i) {
					write(ofp, ";%s", bcf_hdr_int2id(header, BCF_DT_ID, rec->d.flt[i]));
				}
			}
			write(ofp, "\tNS=%d;AN=%d", ns, an);
			write(ofp, ";AC=%d", ac[1]);
			for (int i = 2; i < rec->n_allele; ++i) {
				write(ofp, ",%d", ac[i]);
			}
			write(ofp, ";AF=%g", an > 0 ? (ac[1] / (double)an) : 0.0);
			for (int i = 2; i < rec->n_allele; ++i) {
				write(ofp, ",%g", an > 0 ? (ac[i] / (double)an) : 0.0);
			}
			write(ofp, ";Het=%d", het[1]);
			for (int i = 2; i < rec->n_allele; ++i) {
				write(ofp, ",%d", het[i]);
			}
			write(ofp, ";Hom=%d", hom[1]);
			for (int i = 2; i < rec->n_allele; ++i) {
				write(ofp, ",%d", hom[i]);
			}
			if (!label.empty()) {
				write(ofp, ";VRT=%d\tNA:FRQ\t", vrt);
				write(ofp, "%d:%g", an, an > 0 ? (ac[1] / (double)an) : 0.0);
				for (int i = 2; i < rec->n_allele; ++i) {
					write(ofp, ",%g", an > 0 ? (ac[i] / (double)an) : 0.0);
				}
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
