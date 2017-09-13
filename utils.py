from collections import OrderedDict
from operator import itemgetter
import traceback

AF_BUCKETS = [0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1]
METRICS = {
    'BaseQRankSum':{},
    'ClippingRankSum':{},
    'DP':{'name':'Total Depth'},
    'FS':{},
    'InbreedingCoeff':{},
    'MQ':{'name':'Mapping Quality'},
    'MQRankSum':{},
    'QD':{},
    'ReadPosRankSum':{},
    'VQSLOD':{},
    'SVM':{'name':'SVM Score'},
    'FIBC_P':{'name':'In-Breeding Coefficient'},
    'FIBC_I':{'name':'In-Breeding Coefficient (pop-adjusted)'},
    'HWE_SLP_P':{'name':'HWE p-value'},
    'HWE_SLP_I':{'name':'HWE p-value (pop-adjusted)'},
    'ABE':{'name':'Expected Allele Balance'},
    'ABZ':{'name':'Allele Balance Z-score'},
    'BQZ':{'name':'BaseQual-Allele correlation'},
    'CYZ':{'name':'Cycle-Allele correlation'},
    'STZ':{'name':'Strand-Allele correlation'},
    'IOR':{'name':'Inflated Rate of Observing other alleles (log10)'},
    'NM0':{'name':'Avg num mismatches in reads with ref alleles'},
    'NM1':{'name':'Avg num mismatches in reads with alt alleles'},
    'NMZ':{'name':'Mismatches/read-Allele correlation'},
}
for k,v in METRICS.items(): v.setdefault('name',k)


class Consequence(object):
    # This is a slightly modified version of VEP's recommendations - see http://useast.ensembl.org/info/genome/variation/predicted_data.html#consequences
    # The ordering of the LoF variants is from snpEff's recommendations - see http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf
    # To find all variants that are used, run:
    # mongo --eval 'db.variants.distinct("vep_annotations.Consequence").forEach(printjson)' topmed | tr -d '",' | tr "&" "\n" | sort -u
    _lof_csqs = [
        "transcript_ablation",
        "frameshift_variant",
        "stop_gained",
        "stop_lost",
        "start_lost",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "transcript_amplification",
    ]
    _missense_csqs = [
        "inframe_insertion",
        "inframe_deletion",
        "missense_variant",
        "protein_altering_variant",
    ]
    _synonymous_csqs = [
        "splice_region_variant",
        "incomplete_terminal_codon_variant",
        "stop_retained_variant",
        "synonymous_variant",
    ]
    _other_csqs = [
        "coding_sequence_variant",
        "mature_miRNA_variant",
        "5_prime_UTR_variant",
        "3_prime_UTR_variant",
        "non_coding_transcript_exon_variant",
        "intron_variant",
        "NMD_transcript_variant",
        "non_coding_transcript_variant",
        "upstream_gene_variant",
        "downstream_gene_variant",
        "TFBS_ablation",
        "TFBS_amplification",
        "TF_binding_site_variant",
        "regulatory_region_ablation",
        "regulatory_region_amplification",
        "feature_elongation",
        "regulatory_region_variant",
        "feature_truncation",
        "intergenic_variant",
    ]
    csqs = _lof_csqs + _missense_csqs + _synonymous_csqs + _other_csqs
    assert len(csqs) == len(set(csqs)) # No dupes!
    csqidxs = {csq:i for i,csq in enumerate(csqs)}

class Xpos:
    CHROMOSOME_STRINGS = ['chr%s' % x for x in range(1, 22+1)] + ['chrX', 'chrY', 'chrM']
    CHROMOSOME_STRING_TO_NUMBER = {chrom: idx+1 for idx,chrom in enumerate(CHROMOSOME_STRINGS) }
    CHROMOSOME_NUMBER_TO_STRING = {chrom_num: chrom for chrom,chrom_num in CHROMOSOME_STRING_TO_NUMBER.items()}
    @staticmethod
    def from_chrom_pos(chrom, pos):
        if not chrom.startswith('chr'):
            chrom = 'chr{}'.format(chrom)
        return Xpos.CHROMOSOME_STRING_TO_NUMBER[chrom] * int(1e9) + pos
    @staticmethod
    def to_chrom_pos(xpos):
        pos = xpos % int(1e9)
        chrom = Xpos.CHROMOSOME_NUMBER_TO_STRING[int(xpos) / int(1e9)]
        return (chrom, pos)
    @staticmethod
    def to_pos(xpos):
        return xpos % int(1e9)

class ConsequenceDrilldown(object):
    @staticmethod
    def from_variant(variant):
        """
        Returns something like {"frameshift": {"ENSG00001234": [{"SYMBOL": "APOL1", "Gene": "ENSG00001234", "Feature": "ENST00002345", ...}]}}
        """
        if 'vep_annotations' not in variant:
            return {}
        consequences_drilldown = OrderedDict()
        for annotation in variant['vep_annotations']:
            consequences_drilldown.setdefault(Consequence.csqs[annotation['worst_csqidx']], {}).setdefault(annotation['Gene'], []).append(annotation)
        # Sort the consequences
        for csq in consequences_drilldown:
            for gene in consequences_drilldown[csq]:
                consequences_drilldown[csq][gene] = sorted(consequences_drilldown[csq][gene], key=lambda ann: (ann.get('HGVS'), ann.get('Feature')))
        return consequences_drilldown

    @staticmethod
    def split_into_two_columns(consequences):
            '''
            Try to make two columns of similar height, but with the first a little taller.
            Returns the names of the consequences (ie, the keys), but not the values (because that'd be a pain to use).
            '''
            if len(consequences) == 0:
                return ([], [])
            elif len(consequences) == 1:
                return (consequences.keys(), [])
            consequence_heights = [0]
            for annotations in consequences.values()[0].values():
                consequence_heights[0] += len(annotations) # The number of annotations in this gene (because all are shown in the first consequence)
                # TODO: check for the other things displayed in variant_details.html
            for csq in consequences.values()[1:]:
                consequence_heights.append(len(csq)) # The number of genes in this consequence (because annotations are collapsed in these consequences)
            index = ConsequenceDrilldown._get_midpoint_index(consequence_heights)
            return (consequences.keys()[:index],
                    consequences.keys()[index:])

    @staticmethod
    def _get_midpoint_index(lst):
        '''
        for test_lst in [[1], [1,2,3], [3,1,1], [3,1,1,1], [3,1,1,1,1]]:
            index = get_midpoint_index(test_lst)
            assert 0 < index <= len(test_lst)
            assert sum(test_lst[:index]) >= sum(test_lst[index:])
            assert sum(test_lst[:index-1]) < sum(test_lst[index-1:])
        '''
        half = sum(lst) / 2.0
        acc = 0
        for index, num in enumerate(lst):
            if acc >= half:
                return index
            acc += num
        return len(lst)

    @staticmethod
    def get_top_gene_and_HGVSs(consequences_drilldown):
        """Returns something like ("APOL1", ["Gly70Ter", "Gly88Ter"])"""
        if not consequences_drilldown:
            return None, []
        gene_drilldowns_for_top_csq = consequences_drilldown.values()[0]
        if len(gene_drilldowns_for_top_csq) != 1: # we need exactly one gene
            return None, []
        annotation_drilldowns_for_top_csq = gene_drilldowns_for_top_csq.values()[0]
        gene_symbol_for_top_csq = annotation_drilldowns_for_top_csq[0].get('SYMBOL') or gene_drilldowns_for_top_csq.keys()[0]
        HGVSs_for_top_csq = sorted({ann['HGVS'] for ann in annotation_drilldowns_for_top_csq if ann.get('HGVS')})
        return gene_symbol_for_top_csq, sorted(HGVSs_for_top_csq)

class defaultdict_that_passes_key_to_default_factory(dict):
    "A class like collections.defaultdict, but where the default_factory takes the missing key as an argument."
    def __init__(self, default_factory):
        self._default_factory = default_factory
        super(defaultdict_that_passes_key_to_default_factory, self).__init__()
    def __missing__(self, key):
        value = self[key] = self._default_factory(key)
        return value

def indent_pprint(obj):
    import pprint
    print '\n'.join('####'+line for line in pprint.pformat(obj).split('\n'))
def mkdict(*dicts, **ret):
    for d in dicts: ret.update({k:True for k in d} if isinstance(d, (set,list)) else d)
    return ret

def clamp(num, min_value, max_value):
    return max(min_value, min(max_value, num))

def histogram_from_counter(counter, num_bins=10, bin_range=None):
    from math import floor
    if bin_range is None:
        bin_range = (min(counter.iterkeys()), max(counter.iterkeys()))
    bin_width = float(bin_range[1] - bin_range[0]) / num_bins
    if bin_width == 0:
        only_key = counter.keys()[0]
        print 'Warning: metric always had the value {}'.format(counter.keys())
        return {'left_edges': [only_key-1, only_key, only_key+1], 'mids': [only_key-1, only_key, only_key+1], 'counts': [0, counter.values()[0], 0]}
    bin_left_edges = [bin_range[0] + bin_width * i for i in range(num_bins)]
    bin_counts = [0]*num_bins
    for key, count in counter.iteritems():
        bin_i = (key - bin_range[0]) / bin_width
        try:
            bin_i = int(floor(bin_i))
        except:
            print 'error on', bin_i, key, bin_range[0], bin_range[1], bin_width
            raise
        bin_i = clamp(bin_i, min_value=0, max_value=num_bins-1)
        bin_counts[bin_i] += count
    bin_mids = [left_edge + bin_width/2.0 for left_edge in bin_left_edges]
    return {'left_edges': bin_left_edges, 'mids': bin_mids, 'counts': bin_counts}
