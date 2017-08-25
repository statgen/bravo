import re
from utils import * # TODO: explicitly list
import itertools
import pysam
import json
import time
import pymongo

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
    if variant is None: return None
    if variant['rsids'] == []:
        variant['rsids'] = list('rs{}'.format(r['rsid']) for r in db.dbsnp.find({'xpos': xpos}))
        if variant['rsids']:
            print("apparently the variant [xpos={!r}, ref={!r}, alt={!r}] didn't have any rsids but found some in db.dbsnp")
    variant['genes'] = [gene for gene in variant['genes'] if gene != '']
    return variant

def get_variant_by_variant_id(db, variant_id, default_to_boring_variant=False):
    try:
        chrom, pos, ref, alt = variant_id.split('-')
        pos = int(pos)
        xpos = get_xpos(chrom, pos)
    except:
        return None
    v = get_variant(db, xpos, ref, alt)
    if v is not None:
        return v
    elif default_to_boring_variant:
        return {
            'chrom': chrom,
            'pos': pos,
            'xpos': xpos,
            'ref': ref,
            'alt': alt,
        }
    else:
        return None


def get_variants_by_rsid(db, rsid):
    if not rsid.startswith('rs'):
        return None
    try:
        int(rsid.lstrip('rs'))
    except Exception, e:
        return None
    variants = list(db.variants.find({'rsids': rsid}, projection={'_id': False}))
    for variant in variants:
        add_consequence_to_variant(variant)
        remove_some_extraneous_information(variant)
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
            for variant in variants:
                add_consequence_to_variant(variant)
                remove_some_extraneous_information(variant)
            return variants
    return []

def get_coverage_for_bases(coverages, xstart, xstop=None):
    """
    Get the coverage for the list of bases given by xstart->xstop, inclusive
    Returns list of coverage dicts sorted by pos
    xstop can be None if just one base, but you'll still get back a list
    """
    if xstop is None:
        xstop = xstart
    start_time = time.time()
    coverages_json = coverages.getCoverageX(xstart, xstop)
    print 'tabix\'ed %s base(s) from %s-%s in %s sec' % (len(coverages_json), xstart, xstop, time.time() - start_time)
    return coverages_json


def get_awesomebar_suggestions(autocomplete_strings, query):
    """
    This generates autocomplete suggestions when user
    query is the string that user types
    If it is the prefix for a gene, return list of gene names
    """
    regex = re.compile('^' + re.escape(query), re.IGNORECASE)
    results = (r for r in autocomplete_strings if regex.match(r))
    results = itertools.islice(results, 0, 20)
    return list(results)


# 1:1-1000
_regex_pattern_chr = r'^(?:CHR)?(\d+|X|Y|M|MT)'
_regex_pattern_chr_pos = _regex_pattern_chr + r'\s*[-:/]\s*([\d,]+)'
_regex_pattern_chr_start_end = _regex_pattern_chr_pos + r'\s*[-:/]\s*([\d,]+)'
_regex_pattern_chr_pos_ref_alt = _regex_pattern_chr_pos + r'\s*[-:/]\s*([ATCG]+)\s*[-:/]\s*([ATCG]+)'

_regex_chr = re.compile(_regex_pattern_chr+'$')
_regex_chr_pos = re.compile(_regex_pattern_chr_pos+'$')
_regex_chr_start_end = re.compile(_regex_pattern_chr_start_end+'$')
_regex_chr_pos_ref_alt = re.compile(_regex_pattern_chr_pos_ref_alt+'$')


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

    # Region
    match = _regex_chr.match(query) or _regex_chr_pos.match(query) or _regex_chr_start_end.match(query) or _regex_chr_pos_ref_alt.match(query)
    if match is not None:
        num_groups = len([g for g in match.groups() if g is not None])
        chrom = match.groups()[0]
        if num_groups == 1:
            return 'region', '{}'.format(chrom)
        pos = int(match.groups()[1].replace(',',''))
        if num_groups == 2:
            return 'region', '{}-{}-{}'.format(chrom, pos, pos)
        if num_groups == 3:
            end = int(match.groups()[2].replace(',',''))
            return 'region', '{}-{}-{}'.format(chrom, pos, end)
        return 'variant', '{}-{}-{}-{}'.format(chrom, pos, match.groups()[2], match.groups()[3])

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
    for variant in variants:
        add_consequence_to_variant(variant)
        remove_extraneous_information(variant)
    return list(variants)


def get_metrics(db, variant):
    if 'allele_count' not in variant or variant['allele_num'] == 0:
        return None
    metrics = {}
    for metric in METRICS:
        metrics[metric] = db.metrics.find_one({'metric': metric}, projection={'_id': False})

    # Rename
    for old_name, new_name in [("DP", "Total Depth"), ("MQ", "Mapping Quality")]:
        if old_name in metrics and old_name in variant['quality_metrics']:
            metrics[new_name] = metrics[old_name]
            metrics[new_name]['metric'] = new_name
            del metrics[old_name]
            variant['quality_metrics'][new_name] = variant['quality_metrics'][old_name]
            del variant['quality_metrics'][old_name]
        elif old_name in metrics != old_name in variant['quality_metrics']:
            print('ERROR: problems with {!r} in metrics for variant {!r}'.format(old_name, variant.get('variant_id')))

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


def remove_some_extraneous_information(variant):
    """Remove information not needed by variant.html or any other page"""
    for key in [
            'xpos',
            'xstart',
            'xstop',
            'vep_annotations',
            'pop_acs',
            'pop_ans',
            'pop_homs',
            'sometimes_missense_or_lof',
    ]:
        variant.pop(key, None)

def remove_extraneous_information(variant):
    """Remove information not needed by gene.html, transcript.html or region.html"""
    remove_some_extraneous_information(variant)
    for key in [
            'genotype_depths',
            'genotype_qualities',
            'transcripts',
            'genes',
            'orig_alt_alleles',
            'site_quality',
            'quality_metrics',
    ]:
        variant.pop(key, None)

def get_variants_in_gene(db, gene_id):
    for variant in db.variants.find({'genes': gene_id}, projection={'_id': False}):
        variant['vep_annotations'] = [x for x in variant['vep_annotations'] if x['Gene'] == gene_id]
        add_consequence_to_variant(variant)
        yield variant

def get_most_important_variants_in_gene(db, gene_id):
    variants = []
    for variant in db.variants.find({'genes': gene_id, 'sometimes_missense_or_lof':1}, projection={'_id': False}):
        variant['vep_annotations'] = [x for x in variant['vep_annotations'] if x['Gene'] == gene_id]
        add_consequence_to_variant(variant)
        if variant['category'] in ['lof_variant', 'missense_variant']:
            remove_extraneous_information(variant)
            variants.append(variant)
    return variants

def get_num_variants_in_gene(db, gene_id):
    return db.variants.find({'genes': gene_id}, projection={'_id': False}).count()


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

def get_num_variants_in_transcript(db, transcript_id):
    return db.variants.find({'transcripts': transcript_id}, projection={'_id': False}).count()


def get_exons_in_transcript(db, transcript_id):
    return sorted(list(db.exons.find({'transcript_id': transcript_id, 'feature_type': { "$in": ['CDS', 'UTR', 'exon'] }}, projection={'_id': False})), key=lambda k: k['start'])

def get_exons_in_gene(db, gene_id):
    """
    Returns the "exons", sorted by position.
    """
    return sorted(list(db.exons.find({'gene_id': gene_id, 'feature_type': { "$in": ['CDS', 'UTR', 'exon'] }}, projection={'_id': False})), key=lambda k: k['start'])



def get_variants_for_table(db, chrom, start_pos, end_pos, columns_to_return, order, filter_info, skip, length):
    # 1. match what the user asked for - using [chrom, start_pos, end_pos, filter_info]
    # 2. project to just keys for sorting, sort, and get `_id`s - using [order]
    # 3. get `n_filtered` and `length`-many `_id`s - using [skip, length]
    # 4. look up those `_id`s and project to keys we'll need to make the final result - using [columns_to_return]
    # 5. do modifications ("annotations") in python (LATER: just put [HGVS csq] in mongo) - using [columns_to_return]
    # 6. project to just fields that we return - using [columns_to_return]
    st = time.time()

    mongo_match = [{'xpos': {'$gte': Xpos.from_chrom_pos(chrom, start_pos), '$lte': Xpos.from_chrom_pos(chrom, end_pos)}}]
    if isinstance(filter_info.get('pos_ge',None),int): mongo_match.append({'xpos': {'$gte': Xpos.from_chrom_pos(chrom, filter_info['pos_ge'])}})
    if isinstance(filter_info.get('pos_le',None),int): mongo_match.append({'xpos': {'$lte': Xpos.from_chrom_pos(chrom, filter_info['pos_le'])}})
    if filter_info.get('filter_value',None) is not None:
        if filter_info['filter_value'] == 'PASS': mongo_match.append({'filter': 'PASS'})
        elif filter_info['filter_value'] == 'not PASS': mongo_match.append({'filter': {'$ne': 'PASS'}})
    if isinstance(filter_info.get('maf_ge',None),(float,int)):
        assert 0 <= filter_info['maf_ge'] <= 0.5
        if filter_info['maf_ge'] > 0: mongo_match.append({'$and': [{'allele_freq': {'$gte': filter_info['maf_ge']}},{'allele_freq': {'$lte': 1-filter_info['maf_ge']}}]})
    if isinstance(filter_info.get('maf_le',None),(float,int)):
        assert 0 <= filter_info['maf_le'] <= 0.5
        if filter_info['maf_le'] < 0.5: mongo_match.append({'$or': [{'allele_freq': {'$lte': filter_info['maf_le']}},{'allele_freq': {'$gte': 1-filter_info['maf_le']}}]})

    cols = {
        # after pre-processing, these will look like:
        # <name>: {'sort': {'project': <projection>, 'sort_key': <key>}, 'return': {'project': <projection>, 'annotate': [<annotator>, ...], 'to_client': [<key>, ...]}}
        'allele': {'return': ['rsids', 'ref', 'alt']},
        'pos': {'sort': 'xpos'},
        'hgvs': {'return': {'project': {'vep_annotations':1}, 'annotate': ['add_csq'], 'to_client':['HGVS']}},
        'csq': {'return': {'project': {'vep_annotations':1}, 'annotate': ['add_csq'], 'to_client':['major_consequence']}},
        'filter': {},
        'allele_count': {'sort': True},
        'allele_num': {'sort': True},
        'het': {'sort': {'project': {'het': {'$subtract':['$allele_count',{'$multiply':[2,'$hom_count']}]}}, 'sort_key': 'het'},
                'return': {'project': {'het': {'$subtract':['$allele_count',{'$multiply':[2,'$hom_count']}]}}, 'to_client': ['het']}},
        'hom_count': {'sort': True},
        'allele_freq': {'sort': True},
        'cadd_phred': {'sort': True},
    }
    annotators = {'add_csq': add_consequence_to_variant}
    for name, col in cols.items():
        try:
            if 'sort' not in col: col['sort'] = False
            if col['sort'] == True: col['sort'] = name
            if isinstance(col['sort'], str): col['sort'] = {'project': {col['sort']:1}, 'sort_key':col['sort']}
            assert col['sort'] == False or isinstance(col['sort']['project'], dict) and isinstance(col['sort']['sort_key'], str)
            if 'return' not in col: col['return'] = [name]
            if isinstance(col['return'], list): col['return'] = {'project': {k:1 for k in col['return']}, 'to_client': col['return']}
            assert isinstance(col['return']['project'], dict) and isinstance(col['return']['to_client'], list)
            assert all(ann in annotators for ann in col['return'].get('annotate', []))
        except:
            print('COL = ', col)
            raise

    mongo_projection_before_sort = {}
    mongo_sort = OrderedDict()
    for order_item in order:
        direction = {'asc': pymongo.ASCENDING, 'desc':pymongo.DESCENDING}[order_item['dir']]
        colidx = order_item['column']; colname = columns_to_return[colidx]['name']; col = cols[colname]
        mongo_projection_before_sort.update(col['sort']['project'])
        mongo_sort[col['sort']['sort_key']] = direction

    mongo_projection = mkdict(*[cols[ctr['name']]['return']['project'] for ctr in columns_to_return])
    annotators_to_run = [annotators[ann] for ann in set(itertools.chain.from_iterable(cols[ctr['name']]['return'].get('annotate',[]) for ctr in columns_to_return))]
    keys_to_return = set(itertools.chain.from_iterable(cols[ctr['name']]['return']['to_client'] for ctr in columns_to_return))

    v_ids_curs = db.variants.aggregate([
        {'$match': {'$and': mongo_match}},
        {'$project': mongo_projection_before_sort},
        {'$sort': mongo_sort},
        {'$project': {'_id': 1}},
        {'$group': {'_id':0, 'count':{'$sum':1}, 'results':{'$push':'$$ROOT'}}},
        {'$project': {'_id':0, 'count':1, 'ids':{'$slice':['$results',skip,length]}}},
    ])
    print '## {:0.3f} sec:'.format(time.time()-st), 'cursor created'; st = time.time()
    v_ids_result = list(v_ids_curs)
    if len(v_ids_result) == 0:
        n_filtered, variants = 0, []
    else:
        assert len(v_ids_result) == 1
        n_filtered = v_ids_result[0]['count']
        print '## {:0.3f} sec:'.format(time.time()-st), 'n_filtered={}'.format(n_filtered); st = time.time()
        v_ids = [v['_id'] for v in v_ids_result[0]['ids']]
        variants = [next(db.variants.aggregate([{'$match': {'_id': vid}}, {'$project': mongo_projection}])) for vid in v_ids] # b/c fancy projections require .aggregate()
        print '## {:0.3f} sec:'.format(time.time()-st), 'len(variants)={}'.format(len(variants)); st = time.time()

        for variant in variants:
            for ann in annotators_to_run: ann(variant)
            for key in [key for key in variant if key not in keys_to_return]: del variant[key]
        print '## {:0.3f} sec:'.format(time.time()-st), 'annotated'; st = time.time()

    return {
        'recordsFiltered': n_filtered,
        'recordsTotal': n_filtered,
        'data': variants
    }
