from collections import OrderedDict
from operator import itemgetter

AF_BUCKETS = [0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1]
METRICS = [
    'BaseQRankSum',
    'ClippingRankSum',
    'DP',
    'FS',
    'InbreedingCoeff',
    'MQ',
    'MQRankSum',
    'QD',
    'ReadPosRankSum',
    'VQSLOD'
]


def xpos_to_pos(xpos):
    return int(xpos % 1e9)


def add_consequence_to_variant(variant):
    worst_csq = worst_csq_with_vep(variant['vep_annotations'])
    if worst_csq is None: return
    variant['major_consequence'] = worst_csq['major_consequence']
    variant['HGVSp'] = get_protein_hgvs(worst_csq)
    variant['HGVSc'] = get_transcript_hgvs(worst_csq)
    variant['HGVS'] = get_proper_hgvs(worst_csq)
    variant['CANONICAL'] = worst_csq['CANONICAL']
    if csq_order_dict[variant['major_consequence']] <= csq_order_dict["LOF_THRESHOLD"]:
        variant['category'] = 'lof_variant'
    elif csq_order_dict[variant['major_consequence']] <= csq_order_dict["MISSENSE_THRESHOLD"]:
        variant['category'] = 'missense_variant'
    elif csq_order_dict[variant['major_consequence']] <= csq_order_dict["SYNONYMOUS_THRESHOLD"]:
        variant['category'] = 'synonymous_variant'
    else:
        variant['category'] = 'other_variant'


protein_letters_1to3 = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu',
    'F': 'Phe', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn',
    'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser',
    'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
    'X': 'Ter', '*': 'Ter', 'U': 'Sec'
}


def get_proper_hgvs(annotation):
    # Needs major_consequence
    if annotation['major_consequence'] in ('splice_donor_variant', 'splice_acceptor_variant', 'splice_region_variant'):
        return get_transcript_hgvs(annotation)
    else:
        return get_protein_hgvs(annotation)


def get_transcript_hgvs(annotation):
    return annotation['HGVSc'].split(':')[-1]


def get_protein_hgvs(annotation):
    """
    Takes consequence dictionary, returns proper variant formatting for synonymous variants
    """
    if '%3D' in annotation['HGVSp']: # "%3D" is "="
        try:
            amino_acids = ''.join([protein_letters_1to3[x] for x in annotation['Amino_acids']])
            return "p." + amino_acids + annotation['Protein_position'] + amino_acids
        except Exception, e:
            print 'Could not create HGVS for: %s' % annotation
    return annotation['HGVSp'].split(':')[-1]

# This is a slightly modified version of VEP's recommendations - see http://useast.ensembl.org/info/genome/variation/predicted_data.html#consequences
# The ordering of the LoF variants is from snpEff's recommendations - see http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf
# To find all variants that are used, run:
# mongo --eval 'db.variants.distinct("vep_annotations.Consequence").forEach(printjson)' topmed | tr -d '",' | tr "&" "\n" | sort -u
csq_order = [
    "transcript_ablation",
    "frameshift_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "transcript_amplification",
    "LOF_THRESHOLD",

    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "MISSENSE_THRESHOLD",

    "splice_region_variant",
    "incomplete_terminal_codon_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "SYNONYMOUS_THRESHOLD",

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
assert len(csq_order) == len(set(csq_order)) # No dupes!

csq_order_dict = {csq:i for i,csq in enumerate(csq_order)}
rev_csq_order_dict = dict(enumerate(csq_order))
assert all(csq == rev_csq_order_dict[csq_order_dict[csq]] for csq in csq_order)


def worst_csq_index(csq_list):
    """
    Input list of consequences (e.g. ['frameshift_variant', 'missense_variant'])
    Return index of the worst consequence (In this case, index of 'frameshift_variant', so 4)
    Works well with worst_csq_index('non_coding_exon_variant&nc_transcript_variant'.split('&'))
    """
    # Note: this is the first place that an unknown csq will fail.
    try:
        return min([csq_order_dict[csq] for csq in csq_list])
    except KeyError as exc:
        print("failed on csq_list {!r} with error: {}".format(csq_list, traceback.format_exc()))
        return -9999

def worst_csq_from_list(csq_list):
    """
    Input list of consequences (e.g. ['frameshift_variant', 'missense_variant'])
    Return the worst consequence (In this case, 'frameshift_variant')
    Works well with worst_csq_from_list('non_coding_exon_variant&nc_transcript_variant'.split('&'))
    """
    return rev_csq_order_dict[worst_csq_index(csq_list)]


def worst_csq_from_csq(csq):
    """
    Input possibly &-filled csq string (e.g. 'non_coding_exon_variant&nc_transcript_variant')
    Return the worst consequence (In this case, 'non_coding_exon_variant')
    """
    return rev_csq_order_dict[worst_csq_index(csq.split('&'))]


def order_vep_by_csq(annotation_list):
    """
    Adds "major_consequence" to each annotation.
    Returns them ordered from most deleterious to least.
    """
    for ann in annotation_list:
        ann['major_consequence'] = worst_csq_from_csq(ann['Consequence'])
    return sorted(annotation_list, key=(lambda ann:csq_order_dict[ann['major_consequence']]))


def worst_csq_with_vep(annotation_list):
    """
    Takes list of VEP annotations [{'Consequence': 'frameshift', Feature: 'ENST'}, ...]
    Returns most severe annotation (as full VEP annotation [{'Consequence': 'frameshift', Feature: 'ENST'}])
    Also tacks on "major_consequence" for that annotation (i.e. worst_csq_from_csq)
    """
    if len(annotation_list) == 0:
        return None
    worst = max(annotation_list, key=annotation_severity)
    worst['major_consequence'] = worst_csq_from_csq(worst['Consequence'])
    return worst

def annotation_severity(annotation):
    "Bigger is more important."
    rv = -csq_order_dict[worst_csq_from_csq(annotation['Consequence'])]
    if annotation['CANONICAL'] == 'YES':
        rv += 0.1
    return rv

CHROMOSOMES = ['chr%s' % x for x in range(1, 23)]
CHROMOSOMES.extend(['chrX', 'chrY', 'chrM'])
CHROMOSOME_TO_CODE = { item: i+1 for i, item in enumerate(CHROMOSOMES) }


def get_single_location(chrom, pos):
    """
    Gets a single location from chromosome and position
    chr must be actual chromosme code (chrY) and pos must be integer

    Borrowed from xbrowse
    """
    return CHROMOSOME_TO_CODE[chrom] * int(1e9) + pos


def get_xpos(chrom, pos):
    """
    Borrowed from xbrowse
    """
    if not chrom.startswith('chr'):
        chrom = 'chr{}'.format(chrom)
    return get_single_location(chrom, int(pos))


def get_minimal_representation(pos, ref, alt): 
    """
    Get the minimal representation of a variant, based on the ref + alt alleles in a VCF
    This is used to make sure that multiallelic variants in different datasets, 
    with different combinations of alternate alleles, can always be matched directly. 

    Note that chromosome is ignored here - in xbrowse, we'll probably be dealing with 1D coordinates 
    Args: 
        pos (int): genomic position in a chromosome (1-based)
        ref (str): ref allele string
        alt (str): alt allele string
    Returns: 
        tuple: (pos, ref, alt) of remapped coordinate
    """
    pos = int(pos)
    # If it's a simple SNV, don't remap anything
    if len(ref) == 1 and len(alt) == 1: 
        return pos, ref, alt
    else:
        # strip off identical suffixes
        while(alt[-1] == ref[-1] and min(len(alt),len(ref)) > 1):
            alt = alt[:-1]
            ref = ref[:-1]
        # strip off identical prefixes and increment position
        while(alt[0] == ref[0] and min(len(alt),len(ref)) > 1):
            alt = alt[1:]
            ref = ref[1:]
            pos += 1
        return pos, ref, alt

# Note: I don't know why these are named "MAF".  They are often > 50%.
POP_AFS_1000G = {
    "EAS_MAF": "1000G East Asian",
    "AFR_MAF": "1000G African",
    "EUR_MAF": "1000G European",
    "SAS_MAF": "1000G South Asian",
    "AMR_MAF": "1000G American"
}
def get_pop_afs(variant):
    """
    Convert the nasty output of VEP into a decent dictionary of population AFs.
    """
    if 'vep_annotations' not in variant or len(variant['vep_annotations']) == 0:
        return {}

    pop_strings = {}
    for pop in POP_AFS_1000G:
        values = [ann[pop] for ann in variant['vep_annotations']]
        assert all(value == values[0] for value in values)
        pop_strings[pop] = values[0]

    if all(pop_string == '' for pop_string in pop_strings.values()):
        return {}

    pop_acs = {}
    for pop_key, pop_name in POP_AFS_1000G.items():
        d = {}
        for alt_af in pop_strings[pop_key].split('&'):
            k, v = alt_af.split(':')
            assert all(letter in 'ACTG-' for letter in k)
            d[k] = float(v)
        pop_acs[pop_name] = d.get(variant['alt'])
        if pop_acs[pop_name] is None:
            print('WARNING: pop_af dictionary {!r} is missing alt allele {!r} for population {!r} for variant {}'.format(d, variant['alt'], pop_key, variant['variant_id']))
            return {}
    return pop_acs

def get_consequences_drilldown_for_variant(variant):
    """
    Adds 'major_consequence' to each annotation.
    Returns something like {"frameshift": {"ENSG00001234": [{"SYMBOL": "APOL1", "Gene": "ENSG00001234", "Feature": "ENST00002345", ...}]}}
    """
    if 'vep_annotations' not in variant:
        return {}
    variant['vep_annotations'] = order_vep_by_csq(variant['vep_annotations'])  # Adds major_consequence
    consequences = OrderedDict()
    for annotation in variant['vep_annotations']:
        annotation['HGVS'] = get_proper_hgvs(annotation)
        consequences.setdefault(annotation['major_consequence'], {}).setdefault(annotation['Gene'], []).append(annotation)
    # Sort the consequences
    for csq in consequences:
        for gene in consequences[csq]:
            consequences[csq][gene] = sorted(consequences[csq][gene], key=lambda ann: (ann.get('HGVS'), ann.get('Feature')))
    return consequences

def get_top_gene_and_top_hgvss_for_consequences_drilldown(consequences):
    """Returns something like ("APOL1", ["Gly70Ter", "Gly88Ter"])"""
    if not consequences:
        return None, []
    top_csq = consequences.values()[0]
    if len(top_csq) != 1: # we need exactly one gene
        return None, []
    annotations_for_top_csq = top_csq.values()[0]
    gene_for_top_csq = annotations_for_top_csq[0].get('SYMBOL') or top_csq.keys()[0]
    top_HGVSs = sorted({ann['HGVS'].lstrip('p.') for ann in annotations_for_top_csq if ann.get('HGVS')})
    return gene_for_top_csq, top_HGVSs
