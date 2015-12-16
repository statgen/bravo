import re
from utils import *
import itertools

SEARCH_LIMIT = 10000


def get_gene(db, gene_id):
    return db.genes.find_one({'gene_id': gene_id}, projection={'_id': False})


def get_gene_by_name(db, gene_name):
    # try gene_name field first
    gene = db.genes.find_one({'gene_name': gene_name}, projection={'_id': False})
    if gene:
        return gene
    # if not, try gene['other_names']
    return db.genes.find_one({'other_names': gene_name}, projection={'_id': False})


def get_transcript(db, transcript_id):
    transcript = db.transcripts.find_one({'transcript_id': transcript_id}, projection={'_id': False})
    if not transcript:
        return None
    transcript['exons'] = get_exons_in_transcript(db, transcript_id)
    return transcript


def get_variant(db, xpos, ref, alt):
    variant = db.variants.find_one({'xpos': xpos, 'ref': ref, 'alt': alt}, projection={'_id': False})
    if variant is not None and variant['rsids'] == []:
        variant['rsids'] = list('rs{}'.format(r['rsid']) for r in db.dbsnp.find({'xpos': xpos}))
        if variant['rsids']:
            print("apparently the variant [xpos={!r}, ref={!r}, alt={!r}] didn't have any rsids but found some in db.dbsnp")
    return variant


def get_variants_by_rsid(db, rsid):
    if not rsid.startswith('rs'):
        return None
    try:
        int(rsid.lstrip('rs'))
    except Exception, e:
        return None
    variants = list(db.variants.find({'rsids': rsid}, projection={'_id': False}))
    add_consequence_to_variants(variants)
    return variants


def get_variants_from_dbsnp(db, rsid):
    if not rsid.startswith('rs'):
        return None
    try:
        rsid = int(rsid.lstrip('rs'))
    except Exception, e:
        return None
    position = db.dbsnp.find_one({'rsid': rsid})
    if position:
        variants = list(db.variants.find({'xpos': {'$lte': position['xpos'], '$gte': position['xpos']}}, projection={'_id': False}))
        if variants:
            add_consequence_to_variants(variants)
            return variants
    return []


def get_coverage_for_bases(db, xstart, xstop=None):
    """
    Get the coverage for the list of bases given by xstart->xstop, inclusive
    Returns list of coverage dicts sorted by pos
    xstop can be None if just one base, but you'll still get back a list
    """
    if xstop is None:
        xstop = xstart
    coverages = {
        doc['xpos']: doc for doc in db.base_coverage.find(
            {'xpos': {'$gte': xstart, '$lte': xstop}},
            projection={'_id': False}
        )
    }
    ret = []
    for i in range(xstart, xstop+1):
        if i in coverages:
            ret.append(coverages[i])
        else:
            ret.append({'xpos': i, 'pos': xpos_to_pos(i)})
    for item in ret:
        item['has_coverage'] = 'mean' in item
        del item['xpos']
    return ret


def get_coverage_for_transcript(db, xstart, xstop=None, num_bins=None):
    """
    :param db:
    :param genomic_coord_to_exon:
    :param xstart:
    :param xstop:
    :param num_bins: An approximate intented number of bins.
    :return:
    """
    coverage_array = get_coverage_for_bases(db, xstart, xstop)
    # only return coverages that have coverage (if that makes any sense?)
    # return coverage_array
    covered = [c for c in coverage_array if c['has_coverage']]
    for c in covered:
        del c['has_coverage']

    if num_bins is not None and xstop is not None:
        bin_length = int((xstop - xstart) / num_bins) + 1
        cur_bin = []
        bins = []
        for base in covered:
            if cur_bin == [] or cur_bin[0]['pos']+bin_length > base['pos']:
                cur_bin.append(base)
            else:
                avg_base = {}
                for key in cur_bin[0]:
                    avg_base[key] = sum(b[key] for b in cur_bin) / len(cur_bin)
                avg_base['start_pos'] = cur_bin[0]['pos']
                avg_base['stop_pos'] = cur_bin[-1]['pos']
                del avg_base['pos']
                bins.append(avg_base)
                cur_bin = [base]
        return bins

    return covered


def get_awesomebar_suggestions(g, query):
    """
    This generates autocomplete suggestions when user
    query is the string that user types
    If it is the prefix for a gene, return list of gene names
    """
    regex = re.compile('^' + re.escape(query), re.IGNORECASE)
    results = (r for r in g.autocomplete_strings if regex.match(r))
    results = itertools.islice(results, 0, 20)
    return list(results)


# 1:1-1000
R1 = re.compile(r'^(\d+|X|Y|M|MT)\s*:\s*(\d+)-(\d+)$')
R2 = re.compile(r'^(\d+|X|Y|M|MT)\s*:\s*(\d+)$')
R3 = re.compile(r'^(\d+|X|Y|M|MT)$')
R4 = re.compile(r'^(\d+|X|Y|M|MT)\s*[-:]\s*(\d+)-([ATCG]+)-([ATCG]+)$')


def get_awesomebar_result(db, query):
    """
    Similar to the above, but this is after a user types enter
    We need to figure out what they meant - could be gene, variant, region

    Return tuple of (datatype, identifier)
    Where datatype is one of 'gene', 'variant', or 'region'
    And identifier is one of:
    - ensembl ID for gene
    - variant ID string for variant (eg. 1-1000-A-T)
    - region ID string for region (eg. 1-1000-2000)

    Follow these steps:
    - if query is an ensembl ID, return it
    - if a gene symbol, return that gene's ensembl ID
    - if an RSID, return that variant's string


    Finally, note that we don't return the whole object here - only it's identifier.
    This could be important for performance later

    """
    query = query.strip()
    print 'Query: %s' % query

    # Variant
    variants = get_variants_by_rsid(db, query.lower())
    if variants:
        if len(variants) == 1:
            return 'variant', variants[0]['variant_id']
        else:
            if query.lower() not in variants[0]['rsids']:
                print('Warning: get_variants_by_rsid(db, "{query_lower!r}") returned ({variants!r}) but {query_lower!r} is not in {variants[0].rsids!r}.'.format(
                    query_lower=query.lower(), variants=variants))
            return 'dbsnp_variant_set', query.lower()
    variant = get_variants_from_dbsnp(db, query.lower())
    if variant:
        return 'variant', variant[0]['variant_id']
    # variant = get_variant(db, )
    # TODO - https://github.com/brettpthomas/exac_browser/issues/14

    gene = get_gene_by_name(db, query)
    if gene:
        return 'gene', gene['gene_id']

    # From here out, all should be uppercase (gene, tx, region, variant_id)
    query = query.upper()
    gene = get_gene_by_name(db, query)
    if gene:
        return 'gene', gene['gene_id']

    # Ensembl formatted queries
    if query.startswith('ENS'):
        # Gene
        gene = get_gene(db, query)
        if gene:
            return 'gene', gene['gene_id']

        # Transcript
        transcript = get_transcript(db, query)
        if transcript:
            return 'transcript', transcript['transcript_id']

    # From here on out, only region queries
    if query.startswith('CHR'):
        query = query.lstrip('CHR')
    # Region
    m = R1.match(query)
    if m:
        if int(m.group(3)) < int(m.group(2)):
            return 'region', 'invalid'
        return 'region', '{}-{}-{}'.format(m.group(1), m.group(2), m.group(3))
    m = R2.match(query)
    if m:
        return 'region', '{}-{}-{}'.format(m.group(1), m.group(2), m.group(2))
    m = R3.match(query)
    if m:
        return 'region', '{}'.format(m.group(1))
    m = R4.match(query)
    if m:
        return 'variant', '{}-{}-{}-{}'.format(m.group(1), m.group(2), m.group(3), m.group(4))

    return 'not_found', query


def get_genes_in_region(db, chrom, start, stop):
    """
    Genes that overlap a region
    """
    xstart = get_xpos(chrom, start)
    xstop = get_xpos(chrom, stop)
    genes = db.genes.find({
        'xstart': {'$lte': xstop},
        'xstop': {'$gte': xstart},
    }, projection={'_id': False})
    return list(genes)


def get_variants_in_region(db, chrom, start, stop):
    """
    Variants that overlap a region
    Unclear if this will include CNVs
    """
    xstart = get_xpos(chrom, start)
    xstop = get_xpos(chrom, stop)
    variants = list(db.variants.find({
        'xpos': {'$lte': xstop, '$gte': xstart}
    }, projection={'_id': False}, limit=SEARCH_LIMIT))
    add_consequence_to_variants(variants)
    return list(variants)


def get_metrics(db, variant):
    if 'allele_count' not in variant:
        return None
    metrics = {}
    for metric in METRICS:
        metrics[metric] = db.metrics.find_one({'metric': metric}, projection={'_id': False})

    metric = None
    if variant['allele_count'] == 1:
        metric = 'singleton'
    elif variant['allele_count'] == 2:
        metric = 'doubleton'
    else:
        for af in AF_BUCKETS:
            if float(variant['allele_count'])/variant['allele_num'] < af:
                metric = af
                break
    if metric is not None:
        metrics['Site Quality'] = db.metrics.find_one({'metric': 'binned_%s' % metric}, projection={'_id': False})
    return metrics


def remove_extraneous_information(variant):
    del variant['genotype_depths']
    del variant['genotype_qualities']
    del variant['transcripts']
    del variant['genes']
    del variant['orig_alt_alleles']
    del variant['xpos']
    del variant['xstart']
    del variant['xstop']
    del variant['site_quality']
    del variant['vep_annotations']

def get_variants_in_gene(db, gene_id):
    """
    """
    variants = []
    for variant in db.variants.find({'genes': gene_id}, projection={'_id': False}):
        variant['vep_annotations'] = [x for x in variant['vep_annotations'] if x['Gene'] == gene_id]
        add_consequence_to_variant(variant)
        remove_extraneous_information(variant)
        variants.append(variant)
    return variants

def get_most_important_variants_in_gene(db, gene_id, limit=200):
    # Note: this can almost certainly be heavily optimized.
    lof_variants = []
    missense_variants = []
    other_variants = []
    for variant in db.variants.find({'genes': gene_id}, projection={'_id': False}):
        variant['vep_annotations'] = [x for x in variant['vep_annotations'] if x['Gene'] == gene_id]
        if variant['filter'] != "PASS": continue

        add_consequence_to_variant(variant)
        remove_extraneous_information(variant)

        if variant['category'] == 'lof_variant':
            lof_variants.append(variant)
        elif variant['category'] == 'missense_variant' and len(lof_variants) + len(missense_variants) < limit:
            missense_variants.append(variant)
        elif len(lof_variants) + len(missense_variants) + len(other_variants) < limit:
            other_variants.append(variant)

        if len(lof_variants) == limit:
            return lof_variants
    return lof_variants + missense_variants[:limit - len(lof_variants)] + other_variants[:limit - len(lof_variants) - len(missense_variants)]

def get_num_variants_in_gene(db, gene_id):
    return {
        'total': db.variants.find({'genes': gene_id}, projection={'_id': False}).count(),
        'pass': db.variants.find({'genes': gene_id, 'filter': 'PASS'}, projection={'_id': False}).count()
    }


def get_transcripts_in_gene(db, gene_id):
    """
    """
    return list(db.transcripts.find({'gene_id': gene_id}, projection={'_id': False}))


def get_variants_in_transcript(db, transcript_id):
    """
    """
    variants = []
    for variant in db.variants.find({'transcripts': transcript_id}, projection={'_id': False}):
        variant['vep_annotations'] = [x for x in variant['vep_annotations'] if x['Feature'] == transcript_id]
        add_consequence_to_variant(variant)
        remove_extraneous_information(variant)
        variants.append(variant)
    return variants

def get_most_important_variants_in_transcript(db, transcript_id, limit=200):
    # Note: this can almost certainly be heavily optimized.

    lof_variants = []
    missense_variants = []
    other_variants = []
    for variant in db.variants.find({'transcripts': transcript_id}, projection={'_id': False}):
        variant['vep_annotations'] = [x for x in variant['vep_annotations'] if x['Feature'] == transcript_id]
        if variant['filter'] != "PASS": continue

        add_consequence_to_variant(variant)
        remove_extraneous_information(variant)

        if variant['category'] == 'lof_variant':
            lof_variants.append(variant)
        elif variant['category'] == 'missense_variant' and len(lof_variants) + len(missense_variants) < limit:
            missense_variants.append(variant)
        elif len(lof_variants) + len(missense_variants) + len(other_variants) < limit:
            other_variants.append(variant)

        if len(lof_variants) == limit:
            return lof_variants
    return lof_variants + missense_variants[:limit - len(lof_variants)] + other_variants[:limit - len(lof_variants) - len(missense_variants)]


def get_exons_in_transcript(db, transcript_id):
    # return sorted(
    #     [x for x in
    #      db.exons.find({'transcript_id': transcript_id}, projection={'_id': False})
    #      if x['feature_type'] != 'exon'],
    #     key=lambda k: k['start'])
    return sorted(list(db.exons.find({'transcript_id': transcript_id, 'feature_type': { "$in": ['CDS', 'UTR', 'exon'] }}, projection={'_id': False})), key=lambda k: k['start'])
