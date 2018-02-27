"""
Utils for reading flat files that are loaded into database
"""
import re
import traceback
import itertools
import boltons.iterutils
from utils import *

#Peter: this is no longer used.
# POPS = {
#     'AFR': 'African',
#     'AMR': 'Latino',
#     'EAS': 'East Asian',
#     'FIN': 'European (Finnish)',
#     'NFE': 'European (Non-Finnish)',
#     'SAS': 'South Asian',
#     'OTH': 'Other'
# }


def get_base_coverage_from_file(base_coverage_file):
    """
    Read a base coverage file and return iter of dicts that look like:
    {
        'xpos': 1e9+1,
        'mean': 0.0,
        'median': 0.0,
        '1': 0.0,
        '5': 0.0,
        '10': 0.0,
        '15': 0.0,
        '20': 0.0,
        '25': 0.0,
        '30': 0.0,
        '50': 0.0,
        '100': 0.0,
    }
    """

    float_header_fields = ['mean', 'median', '1', '5', '10', '15', '20', '25', '30', '50', '100']
    for line in base_coverage_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')
        d = {
            'xpos': Xpos.from_chrom_pos(fields[0], int(fields[1])),
            'pos': int(fields[1]),
        }
        for i, k in enumerate(float_header_fields):
            d[k] = float(fields[i+2])
        yield d

def get_variants_from_sites_vcf_without_annotation(sites_vcf):
    for line in sites_vcf:
        try:
            line = line.strip('\n')
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            info_field = dict(x.split('=', 1) if '=' in x else (x, x) for x in re.split(';(?=\w)', fields[7]))
            for i, alt_allele in enumerate(fields[4].split(',')):
                variant = {}
                variant['chrom'] = fields[0][3:] if fields[0].startswith('chr') else fields[0]
                variant['pos'], variant['ref'], variant['alt'] = get_minimal_representation(fields[1], fields[3], alt_allele)
                variant['xpos'] = Xpos.from_chrom_pos(variant['chrom'], variant['pos'])
                variant['xstop'] = variant['xpos'] + len(variant['alt']) - len(variant['ref'])
                variant['variant_id'] = '{}-{}-{}-{}'.format(variant['chrom'], variant['pos'], variant['ref'], variant['alt'])
                variant['site_quality'] = float(fields[5])
                variant['filter'] = fields[6]
                variant['allele_count'] = int(info_field['AC'].split(',')[i])
                if variant['allele_count'] == 0: continue
                variant['allele_num'] = int(info_field['AN'])
                assert variant['allele_num'] != 0, variant
                variant['allele_freq'] = float(info_field['AF'].split(',')[i])
                assert variant['allele_freq'] != 0, variant
                variant['hom_count'] = int(info_field['Hom'].split(',')[i])
                yield variant
        except:
            print("Error while parsing vcf line: " + line)
            traceback.print_exc()
            raise


def get_variants_from_sites_vcf_only_percentiles(sites_vcf):
    for line in sites_vcf:
        try:
            line = line.strip('\n')
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            percentiles = {}
            for info_entry in re.split(';(?=\w)', fields[7]):
                info_entry = info_entry.split('=')
                if len(info_entry) == 2 and info_entry[0].endswith('_PCTL'):
                    value = map(float, info_entry[1].split(','))
                    if len(value) == 2:
                        key = info_entry[0][:-5]
                        percentiles[key] = value
            if not percentiles:
                continue
            for alt in fields[4].split(','):
                variant = {}
                variant['chrom'] = fields[0][3:] if fields[0].startswith('chr') else fields[0]
                variant['pos'], variant['ref'], variant['alt'] = get_minimal_representation(fields[1], fields[3], alt)
                variant['xpos'] = Xpos.from_chrom_pos(variant['chrom'], variant['pos'])
                variant['percentiles'] = percentiles
                yield variant
        except:
            print("Error while parsing vcf line: " + line)
            traceback.print_exc()
            raise


def get_variants_from_sites_vcf(sites_vcf):
    """
    Parse exac sites VCF file and return iter of variant dicts
    sites_vcf is a file (gzipped), not file path
    """

    # DT: dirty temporary fill-in:
    #dp_mids = map(float, '2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5'.split('|'))
    #gq_mids = map(float, '2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5'.split('|'))
    #hists_all = ['0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0', '0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0']

    vep_field_names = None
    for line in sites_vcf:
        try:
            line = line.strip('\n')
            if line.startswith('##INFO=<ID=CSQ'):
                vep_field_names = line.split('Format: ')[-1].strip('">').split('|')
            if line.startswith('##INFO=<ID=DP_HIST'):
                dp_mids = [float(x) for x in line.split('Mids: ')[-1].strip('">').split('|')]
            if line.startswith('##INFO=<ID=GQ_HIST'):
                gq_mids = [float(x) for x in line.split('Mids: ')[-1].strip('">').split('|')]
            if line.startswith('#'):
                continue

            # It's a variant line

            if vep_field_names is None:
                raise Exception("VEP_field_names is None. Make sure VCF header is present.")
            # This elegant parsing code below is copied from https://github.com/konradjk/loftee
            fields = line.split('\t')
            info_field = dict(x.split('=', 1) if '=' in x else (x, x) for x in re.split(';(?=\w)', fields[7]))

            if 'DP_HIST' not in info_field: continue # This skips weird lines introduced by some kind of bug.  TODO: fix that bug.

            annotation_arrays = [anno_array.split('|') for anno_array in info_field['CSQ'].split(',')] if 'CSQ' in info_field else []
            for anno_array in annotation_arrays: assert len(vep_field_names) == len(anno_array), (vep_field_names, anno_array)
            annotations = [dict(zip(vep_field_names, anno_array)) for anno_array in annotation_arrays]

            # make a different variant for each alt allele
            # we could do the parsing all at once, but multi-allelic variants are rare so who cares
            for i, alt_allele in enumerate(fields[4].split(',')):

                vep_annotations = [ann for ann in annotations if int(ann['ALLELE_NUM']) == i + 1]

                variant = {}
                variant['chrom'] = fields[0][3:] if fields[0].startswith('chr') else fields[0]
                variant['pos'], variant['ref'], variant['alt'] = get_minimal_representation(fields[1], fields[3], alt_allele)

                variant['cadd_raw'] = float(info_field['CADD_RAW'].split(',')[i]) if 'CADD_RAW' in info_field else None
                variant['cadd_phred'] = float(info_field['CADD_PHRED'].split(',')[i]) if 'CADD_PHRED' in info_field else None
                if 'AVGDP' in info_field: variant['avgdp'] = float(info_field['AVGDP'])

                #DT: We will take rsIds from the VEP annotation
                #variant['rsid'] = fields[2]
                if vep_annotations:
                   variant['rsids'] = [rsid for rsid in vep_annotations[0]['Existing_variation'].split('&') if rsid.startswith("rs")]
                else:
                   variant['rsids'] = []

                variant['xpos'] = Xpos.from_chrom_pos(variant['chrom'], variant['pos'])
                variant['xstop'] = variant['xpos'] + len(variant['alt']) - len(variant['ref'])
                variant['variant_id'] = '{}-{}-{}-{}'.format(variant['chrom'], variant['pos'], variant['ref'], variant['alt'])
                variant['site_quality'] = float(fields[5])
                variant['filter'] = fields[6]
                variant['vep_annotations'] = vep_annotations

                # DT: variant['allele_count'] = int(info_field['AC_Adj'].split(',')[i])
                variant['allele_count'] = int(info_field['AC'].split(',')[i])
                if variant['allele_count'] == 0: continue
                # DT: if not variant['allele_count'] and variant['filter'] == 'PASS': variant['filter'] = 'AC_Adj0' # Temporary filter
                # DT: variant['allele_num'] = int(info_field['AN_Adj'])
                variant['allele_num'] = int(info_field['AN'])
                assert variant['allele_num'] != 0, variant
                # DT: if variant['allele_num'] > 0:
                # DT:    variant['allele_freq'] = variant['allele_count']/float(info_field['AN_Adj'])
                # DT: else:
                variant['allele_freq'] = float(info_field['AF'].split(',')[i])
                assert variant['allele_freq'] != 0, variant

                #PJVH: these are always zero, so just use get_pop_afs() below.
                # # DT: variant['pop_acs'] = dict([(POPS[x], int(info_field['AC_%s' % x].split(',')[i])) for x in POPS])
                # variant['pop_acs'] = dict([(POPS[x], 0) for x in POPS])
                # # DT: variant['pop_ans'] = dict([(POPS[x], int(info_field['AN_%s' % x])) for x in POPS])
                # variant['pop_ans'] = dict([(POPS[x], 0) for x in POPS])
                # # DT: variant['pop_homs'] = dict([(POPS[x], int(info_field['Hom_%s' % x].split(',')[i])) for x in POPS])
                # variant['pop_homs'] = dict([(POPS[x], 0) for x in POPS])

                # DT: variant['hom_count'] = sum(variant['pop_homs'].values())
                variant['hom_count'] = int(info_field['Hom'].split(',')[i])
                #if variant['chrom'] in ('X', 'Y'):
                    # DT: variant['pop_hemis'] = dict([(POPS[x], int(info_field['Hemi_%s' % x].split(',')[i])) for x in POPS])
                    #PJVH: variant['pop_hemis'] = dict([(POPS[x], 0) for x in POPS])
                    #PJVH: variant['hemi_count'] = sum(variant['pop_hemis'].values())

                variant['quality_metrics'] = {x: info_field[x] for x in METRICS if x in info_field}

                variant['genes'] = list(set(annotation['Gene'] for annotation in vep_annotations if annotation['Gene']))
                variant['transcripts'] = list(set(annotation['Feature'] for annotation in vep_annotations if annotation['Feature']))

                if 'DP_HIST' in info_field:
                   hists_all = [info_field['DP_HIST'].split(',')[0], info_field['DP_HIST'].split(',')[i+1]]
                   variant['genotype_depths'] = [zip(dp_mids, map(int, x.split('|'))) for x in hists_all]
                if 'GQ_HIST' in info_field:
                   hists_all = [info_field['GQ_HIST'].split(',')[0], info_field['GQ_HIST'].split(',')[i+1]]
                   variant['genotype_qualities'] = [zip(gq_mids, map(int, x.split('|'))) for x in hists_all]
                
                # DT: dirty temporary fill-in:
                # DT: variant['genotype_depths'] = [zip(dp_mids, map(int, x.split('|'))) for x in hists_all]
                # DT: variant['genotype_qualities'] = [zip(gq_mids, map(int, x.split('|'))) for x in hists_all] 

                clean_annotation_consequences_for_variant(variant)
                pop_afs = get_pop_afs(variant)
                if pop_afs: variant['pop_afs'] = pop_afs
                keep_only_needed_annotation_fields(variant)
                yield variant
        except:
            print("Error parsing vcf line: " + line)
            traceback.print_exc()
            raise

def get_minimal_representation(pos, ref, alt):
    """
    Get the minimal representation of a variant, based on the ref + alt alleles in a VCF
    This is used to make sure that multiallelic variants in different datasets,
    with different combinations of alternate alleles, can always be matched directly.
    """
    pos = int(pos)
    # If it's a simple SNV, don't remap anything
    if len(ref) == 1 and len(alt) == 1:
        return (pos, ref, alt)
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
        return (pos, ref, alt)

def clean_annotation_consequences_for_variant(variant):
    '''
    add variant.vep_annotions[*].{HGVS,worst_csqidx}
    sort variant.vep_annotations by severity.
    add variant.worst_csq*.
    '''
    if len(variant['vep_annotations']) == 0:
        raise Exception('why no annos for {!r}?'.format(variant))

    for anno in variant['vep_annotations']:
        anno['CANONICAL'] = (anno['CANONICAL'] == 'YES')
        anno['worst_csqidx'] = _get_worst_csqidx_for_annotation(anno)
        anno['HGVS'] = _get_hgvs(anno)
    variant['vep_annotations'] = sorted(variant['vep_annotations'], key=_annotation_severity, reverse=True)

    # TODO: should we just query on `variant.vep_annotations[0].worst_csqidx`? or make `variant.worst_annotation = variant.vep_annotations[0]`?
    worst_anno = variant['vep_annotations'][0]
    variant['worst_csq_CANONICAL'] = worst_anno['CANONICAL']
    variant['worst_csqidx'] = worst_anno['worst_csqidx']
    variant['worst_csq_HGVS'] = worst_anno['HGVS']

def _get_worst_csqidx_for_annotation(annotation):
    try:
        return min(Consequence.csqidxs[csq] for csq in annotation['Consequence'].split('&'))
    except KeyError:
        raise Exception("failed to get csqidx for {!r} with error: {}".format(annotation['Consequence'], traceback.format_exc()))
def _annotation_severity(annotation):
    "higher is more deleterious"
    rv = -annotation['worst_csqidx']
    if annotation['CANONICAL']: rv += 0.1
    return rv
def _get_hgvs(annotation):
    # ExAC code did fancy things, but this feels okay to me.
    from urllib import unquote
    hgvsp = unquote(annotation['HGVSp']).split(':',1)[-1]
    hgvsc = unquote(annotation['HGVSc']).split(':',1)[-1]
    if hgvsp and '=' not in hgvsp: return hgvsp
    if hgvsc: return hgvsc
    if hgvsp: return hgvsp
    return ''

POP_AFS_1000G = {
    "EAS_AF": "1000G East Asian",
    "AFR_AF": "1000G African",
    "EUR_AF": "1000G European",
    "SAS_AF": "1000G South Asian",
    "AMR_AF": "1000G American"
}
def get_pop_afs(variant):
    """
    Convert the nasty output of VEP into a decent dictionary of population AFs.
    """
    if 'vep_annotations' not in variant or len(variant['vep_annotations']) == 0:
        return {}
    try:
        pop_afs = {}
        for pop in POP_AFS_1000G:
            af_strings = [ann[pop].split('&')[0] for ann in variant['vep_annotations'] if ann['Allele'] == variant['alt'] or ann['Allele'] == '-']
            assert boltons.iterutils.same(af_strings)
            if af_strings and af_strings[0] != '':
                pop_afs[pop] = float(af_strings[0])
        if all(v==0 for v in pop_afs.values()):
            return {}
        return pop_afs
    except:
        print('failed in get_pop_afs() for variant {!r} with error:'.format(variant, traceback.format_exc()))
        raise

def keep_only_needed_annotation_fields(variant):
    pass # keep everything for now.
    # needed_fields = ['CANONICAL', 'Gene', 'SYMBOL', 'worst_csqidx', 'Feature']
    # for anno in variant['vep_annotations']:
    #     for unused_field in [field for field in anno if field not in needed_fields]:
    #         del anno[unused_field]



def get_canonical_transcripts(canonical_transcript_file):
    for line in canonical_transcript_file:
        gene, transcript = line.strip().split()
        yield gene, transcript


def get_omim_associations(omim_file):
    for line in omim_file:
        fields = line.strip().split('\t')
        if len(fields) == 4:
            yield fields
        else:
            yield None


def get_genes_from_gencode_gtf(gtf_file):
    """
    Parse gencode GTF file;
    Returns iter of gene dicts
    """
    for line in gtf_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')

        if fields[2] != 'gene':
            continue

        chrom = fields[0][3:]
        start = int(fields[3]) + 1  # bed files are 0-indexed
        stop = int(fields[4]) + 1
        info = dict(x.strip().split() for x in fields[8].split(';') if x != '')
        info = {k: v.strip('"') for k, v in info.items()}
        gene_id = info['gene_id'].split('.')[0]

        gene = {
            'gene_id': gene_id,
            'gene_name': info['gene_name'],
            'gene_name_upper': info['gene_name'].upper(),
            'chrom': chrom,
            'start': start,
            'stop': stop,
            'strand': fields[6],
            'xstart': Xpos.from_chrom_pos(chrom, start),
            'xstop': Xpos.from_chrom_pos(chrom, stop),
        }
        yield gene


def get_transcripts_from_gencode_gtf(gtf_file):
    """
    Parse gencode GTF file;
    Returns iter of transcript dicts
    """
    for line in gtf_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')

        if fields[2] != 'transcript':
            continue

        chrom = fields[0][3:]
        start = int(fields[3]) + 1  # bed files are 0-indexed
        stop = int(fields[4]) + 1
        info = dict(x.strip().split() for x in fields[8].split(';') if x != '')
        info = {k: v.strip('"') for k, v in info.items()}
        transcript_id = info['transcript_id'].split('.')[0]
        gene_id = info['gene_id'].split('.')[0]

        gene = {
            'transcript_id': transcript_id,
            'gene_id': gene_id,
            'chrom': chrom,
            'start': start,
            'stop': stop,
            'strand': fields[6],
            'xstart': Xpos.from_chrom_pos(chrom, start),
            'xstop': Xpos.from_chrom_pos(chrom, stop),
        }
        yield gene


def get_exons_from_gencode_gtf(gtf_file):
    """
    Parse gencode GTF file;
    Returns iter of transcript dicts
    """
    for line in gtf_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')

        if fields[2] not in ['exon', 'CDS', 'UTR']:
            continue

        chrom = fields[0][3:]
        feature_type = fields[2]
        start = int(fields[3]) + 1  # bed files are 0-indexed
        stop = int(fields[4]) + 1
        info = dict(x.strip().split() for x in fields[8].split(';') if x != '')
        info = {k: v.strip('"') for k, v in info.items()}
        transcript_id = info['transcript_id'].split('.')[0]
        gene_id = info['gene_id'].split('.')[0]

        exon = {
            'feature_type': feature_type,
            'transcript_id': transcript_id,
            'gene_id': gene_id,
            'chrom': chrom,
            'start': start,
            'stop': stop,
            'strand': fields[6],
            'xstart': Xpos.from_chrom_pos(chrom, start),
            'xstop': Xpos.from_chrom_pos(chrom, stop),
        }
        yield exon


def get_dbnsfp_info(dbnsfp_file):
    """
    Parse dbNSFP_gene file;
    Returns iter of transcript dicts
    """
    header = dbnsfp_file.next().split('\t')
    fields = dict(zip(header, range(len(header))))
    for line in dbnsfp_file:
        line = line.split('\t')
        other_names = line[fields["Gene_old_names"]].split(';') if line[fields["Gene_old_names"]] != '.' else []
        if line[fields["Gene_other_names"]] != '.':
            other_names.extend(line[fields["Gene_other_names"]].split(';'))
        gene_info = {
            'gene_name': line[fields["Gene_name"]],
            'ensembl_gene': line[fields["Ensembl_gene"]],
            'gene_full_name': line[fields["Gene_full_name"]],
            'gene_other_names': other_names
        } 
        yield gene_info


def get_snp_from_dbsnp_file(dbsnp_file):
    for line in dbsnp_file:
        fields = line.split('\t')
        if len(fields) < 3: continue
        rsid = int(fields[0])
        chrom = fields[1].rstrip('T')
        if chrom == 'PAR': continue
        start = int(fields[2]) + 1
        snp = {
            'xpos': Xpos.from_chrom_pos(chrom, start),
            'rsid': rsid
        }
        yield snp
