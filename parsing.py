"""
Utils for reading flat files that are loaded into database
"""
import itertools
import re
import traceback
from contextlib import closing

import boltons.iterutils
import pysam
from utils import *


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


def get_variants_from_sites_vcf(vcf, chrom, start_bp, end_bp, histograms = True):
    """Reads sites VCF/BCF file and returns iterator over veriant dicts.

    Arguments:
    vcf -- VCF/BCF file name.
    chrom -- chromosome name.
    start_bp -- start position in base-pairs.
    end_bp -- end position in base-pairs.
    histograms -- if True, includes DP and GQ histograms.
    """
    with closing(pysam.VariantFile(vcf)) as ifile:
        vep_meta = ifile.header.info.get('CSQ', None)
        if vep_meta is None:
            raise Exception('Missing CSQ INFO field from VEP (Variant Effect Predictor)')
        vep_field_names = vep_meta.description.split(':', 1)[-1].strip().split('|')
        for x in ['AVGDP', 'AVGDP_R', 'AVGGQ', 'AVGGQ_R']:
            if x not in ifile.header.info:
                raise Exception('Missing {} INFO field.'.format(x))
        if histograms:
            for x in ['DP_HIST', 'DP_HIST_R', 'GQ_HIST', 'GQ_HIST_R']:
                meta = ifile.header.info.get(x, None)
                if x not in ifile.header.info:
                    raise Exception('Missing {} INFO field.'.format(x))
            dp_hist_mids = map(float, ifile.header.info['DP_HIST'].description.split(':', 1)[-1].strip().split('|'))
            dp_hist_r_mids = map(float, ifile.header.info['DP_HIST_R'].description.split(':', 1)[-1].strip().split('|'))
            gq_hist_mids = map(float, ifile.header.info['GQ_HIST'].description.split(':', 1)[-1].strip().split('|'))
            gq_hist_r_mids = map(float, ifile.header.info['GQ_HIST_R'].description.split(':', 1)[-1].strip().split('|'))
        for record in ifile.fetch(chrom, start_bp, end_bp):
            try:
                annotations = dict()
                for annotation in record.info['CSQ']:
                    annotation = annotation.split('|')
                    assert len(vep_field_names) == len(annotation), (vep_field_names, annotation)
                    annotation = dict(zip(vep_field_names, annotation))
                    annotations.setdefault(int(annotation['ALLELE_NUM']), []).append(annotation)
                for i, alt_allele in enumerate(record.alts):
                    variant = {}
                    variant['chrom'] = record.contig[3:] if record.contig.startswith('chr') else record.contig
                    variant['pos'], variant['ref'], variant['alt'] = get_minimal_representation(record.pos, record.ref, alt_allele)
                    variant['xpos'] = Xpos.from_chrom_pos(variant['chrom'], variant['pos'])
                    variant['xstop'] = variant['xpos'] + len(variant['alt']) - len(variant['ref'])
                    variant['variant_id'] = '{}-{}-{}-{}'.format(variant['chrom'], variant['pos'], variant['ref'], variant['alt'])
                    allele_annotations = annotations[i + 1]
                    if allele_annotations:
                        variant['rsids'] = [rsid for rsid in allele_annotations[0]['Existing_variation'].split('&') if rsid.startswith('rs')]
                    else:
                        variant['rsids'] = []
                    variant['site_quality'] = record.qual
                    variant['filter'] = ';'.join(record.filter.keys())
                    variant['allele_count'] = record.info['AC'][i]
                    if variant['allele_count'] == 0:
                        continue
                    variant['allele_num'] = record.info['AN']
                    assert variant['allele_num'] != 0, variant
                    variant['allele_freq'] = record.info['AF'][i]
                    assert variant['allele_freq'] != 0, variant
                    variant['hom_count'] = record.info['Hom'][i]
                    variant['quality_metrics'] = {x: record.info[x] for x in METRICS if x in record.info}
                    variant['genes'] = list(set(annotation['Gene'] for annotation in allele_annotations if annotation['Gene']))
                    variant['transcripts'] = list(set(annotation['Feature'] for annotation in allele_annotations if annotation['Feature']))
                    variant['avgdp'] = record.info['AVGDP']
                    variant['avgdp_alt'] = record.info['AVGDP_R'][i + 1]
                    variant['avggq'] = record.info['AVGGQ']
                    variant['avggq_alt'] = record.info['AVGGQ_R'][i + 1]
                    variant['cadd_raw'] = record.info['CADD_RAW'][i] if 'CADD_RAW' in record.info else None
                    variant['cadd_phred'] = record.info['CADD_PHRED'][i] if 'CADD_PHRED' in record.info else None
                    if histograms:
                        variant['genotype_depths'] = map(lambda x, y: zip(x, map(int, y.split('|'))), [dp_hist_mids, dp_hist_r_mids], [record.info['DP_HIST'], record.info['DP_HIST_R'][i + 1]])
                        variant['genotype_qualities'] = map(lambda x, y: zip(x, map(int, y.split('|'))), [gq_hist_mids, gq_hist_r_mids], [record.info['GQ_HIST'], record.info['GQ_HIST_R'][i + 1]])
                    variant['vep_annotations'] = allele_annotations
                    clean_annotation_consequences_for_variant(variant)
                    pop_afs = get_pop_afs(variant)
                    if pop_afs:
                        variant['pop_afs'] = pop_afs
                    keep_only_needed_annotation_fields(variant)
                    yield variant
            except:
                print("Error parsing VCF/BCF record: " + record.__str__())
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
    header = omim_file.readline().rstrip('\n').split('\t')
    for line in omim_file:
        fields = line.rstrip('\n').split('\t')
        assert len(header) == len(fields)
        fields = dict(zip(header, fields))
        if not fields['MIM gene accession'] or not fields['MIM gene description']:
            continue
        yield fields['Gene stable ID'], fields['Transcript stable ID'], fields['MIM gene accession'], fields['MIM gene description']


def get_regions_from_gencode_gtf(gtf_file, region_types):
    """
    Parse gencode GTF file.
    Returns iter of regions ditcs
    """
    for line in gtf_file:
        if line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        if fields[2] not in region_types:
            continue
        chrom = fields[0][3:]
        start = long(fields[3])
        stop = long(fields[4])
        info = dict(x.strip().split() for x in fields[8].split(';') if x != '')
        region = {
            'chrom': chrom,
            'start': start,
            'stop': stop,
            'strand': fields[6],
            'xstart': Xpos.from_chrom_pos(chrom, start),
            'xstop': Xpos.from_chrom_pos(chrom, stop),
            'gene_id': info['gene_id'].strip('"').split('.')[0],
        }
        if 'gene' in region_types:
            region['gene_name'] = info['gene_name'].strip('"')
        if 'transcript' in region_types:
            region['transcript_id'] = info['transcript_id'].strip('"').split('.')[0] if 'transcript_id' in info else None
        if 'exon' in region_types or 'CDS' in region_types or 'UTR' in region_types:
            if 'transcript_id' not in region:
                region['transcript_id'] = info['transcript_id'].strip('"').split('.')[0] if 'transcript_id' in info else None
            region['feature_type'] = fields[2]
        yield region


def get_genenames(genenames_file):
    """
    Parse file with genes from HGNC.
    Returns iter of gene dicts.
    """
    header = genenames_file.readline().strip('\n').split('\t')
    for line in genenames_file:
        fields = line.rstrip('\n').split('\t')
        assert len(header) == len(fields)
        fields = dict(zip(header, fields))
        if not fields['ensembl_gene_id']:
            continue
        gene = {
            'gene_name': fields['symbol'],
            'ensembl_gene': fields['ensembl_gene_id'],
            'gene_full_name': fields['name'],
            'gene_other_names': fields['alias_symbol'].strip('"').split('|') if fields['alias_symbol'] else []
        }
        if fields['prev_symbol'] and fields['prev_symbol'] not in gene['gene_other_names']:
            for name in fields['prev_symbol'].strip('"').split('|'):
                gene['gene_other_names'].append(name)
        yield gene


def get_snp_from_dbsnp_file(dbsnp_file, chrom, start_bp = None, end_bp = None, histograms = False):
    with closing(pysam.Tabixfile(dbsnp_file, 'r')) as tabix:
        chrom_out = chrom[:-1] if chrom == 'MT' or chrom == 'chrMT' else chrom
        for row in tabix.fetch(chrom, start_bp, end_bp):
            fields = row.split('\t')
            if len(fields) == 3:
                yield { 'xpos': Xpos.from_chrom_pos(chrom_out, int(fields[2]) + 1), 'rsid': int(fields[0]) }
