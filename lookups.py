import itertools
import json
import re
import time

import boltons.iterutils
import pymongo
import pysam
from utils import *  # TODO: explicitly list

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
    return db.transcripts.find_one({'transcript_id': transcript_id}, projection={'_id': False})


def get_variant(db, xpos, ref, alt):
    variant = db.variants.find_one({'xpos': xpos, 'ref': ref, 'alt': alt}, projection={'_id': False})
    if variant is None: return None
    if variant['rsids'] == []:
        variant['rsids'] = list('rs{}'.format(r['rsid']) for r in db.dbsnp.find({'xpos': xpos}))
        if variant['rsids']:
            print("apparently the variant [xpos={!r}, ref={!r}, alt={!r}] didn't have any rsids but found some in db.dbsnp".format(xpos, ref, alt))
    variant['genes'] = [gene for gene in variant['genes'] if gene != '']
    return variant

def get_variant_by_variant_id(db, variant_id, default_to_boring_variant = False):
    try:
        chrom, pos, ref, alt = variant_id.split('-')
        pos = int(pos)
        xpos = Xpos.from_chrom_pos(chrom, pos)
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
    if not rsid.startswith('rs') or not rsid[2:].isdigit():
        return None
    variants = list(db.variants.find({'rsids': rsid}, projection={'_id': False}))
    for variant in variants:
        remove_some_extraneous_information(variant)
    return variants


def get_variants_from_dbsnp(db, rsid):
    if not rsid.startswith('rs') or not rsid[2:].isdigit():
        return None
    position = db.dbsnp.find_one({'rsid': rsid})
    if position:
        variants = list(db.variants.find({'xpos': {'$lte': position['xpos'], '$gte': position['xpos']}}, projection={'_id': False}))
        if variants:
            for variant in variants:
                remove_some_extraneous_information(variant)
            return variants
    return []

#@boltons.cacheutils.cached({})

def get_awesomebar_suggestions(autocomplete_strings, query, db):
    cap = 10
    rs_max_length = 9999999999

    # first look for genes, genes have priority over rsIds (e.g. there is a gene RS1)
    regex = re.compile('^' + re.escape(query), re.IGNORECASE)
    results = (r for r in autocomplete_strings if regex.match(r))
    results = list(itertools.islice(results, 0, cap))

    try:
        if len(results) < cap and query.startswith('rs'): # if query starts with "rs" and there is still place for autocomplete dropdown, look for rsIds.
            rs_numeric = int(query[2:]) if len(query) > 2 else 0
            results.extend('rs{}'.format(x['rsid']) for x  in db.dbsnp.find({ 'rsid': { '$eq': rs_numeric }}, projection = { '_id': False, 'rsid': True }).limit(cap - len(results)))
            step = 10
            while len(results) < cap:
                rs_numeric_lo = rs_numeric * step
                rs_numeric_up = rs_numeric_lo + step - 1
                if rs_numeric_up > rs_max_length:
                    break
                results.extend('rs{}'.format(x['rsid']) for x  in db.dbsnp.find({ 'rsid': { '$gte': rs_numeric_lo, '$lte': rs_numeric_up }}, projection = { '_id': False, 'rsid': True }).limit(cap - len(results)))
                step *= 10
    except ValueError:
        pass
    return results


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
    query = query.strip() # TODO:check if query is not None

    # rsid
    variants = get_variants_by_rsid(db, query.lower())
    if variants:
        if len(variants) == 1:
            return 'variant', {'variant_id': variants[0]['variant_id']}
        else:
            if query.lower() not in variants[0]['rsids']:
                print('Warning: get_variants_by_rsid(db, "{query_lower!r}") returned ({variants!r}) but {query_lower!r} is not in {variants[0].rsids!r}.'.format(
                    query_lower=query.lower(), variants=variants))
            return 'multi_variant_rsid', {'rsid': query.lower()}
    variants = get_variants_from_dbsnp(db, query.lower())
    if variants:
        if len(variants) == 1:
            return 'variant', {'variant_id': variants[0]['variant_id']}
        else:
            return 'multi_variant_rsid', {'rsid': query.lower()}

    # gene symbol
    gene = get_gene_by_name(db, query)
    if gene:
        return 'gene', {'gene_id': gene['gene_id']}

    # From here out, all should be uppercase (gene, tx, region, variant_id)
    query_orig = query
    query = query.upper()

    # uppercase gene symbol
    gene = get_gene_by_name(db, query)
    if gene:
        return 'gene', {'gene_id': gene['gene_id']}

    # ENSG
    if query.startswith('ENSG'):
        gene = get_gene(db, query)
        if gene:
            return 'gene', {'gene_id': gene['gene_id']}

    # ENST
    if query.startswith('ENST'):
        transcript = get_transcript(db, query)
        if transcript:
            return 'transcript', {'transcript_id': transcript['transcript_id']}

    # Region (chrom , chrom-pos , chrom-start-stop) or Variant (chrom-pos-ref-alt)
    match = _regex_chr.match(query) or _regex_chr_pos.match(query) or _regex_chr_start_end.match(query) or _regex_chr_pos_ref_alt.match(query)
    if match is not None:
        num_groups = len([g for g in match.groups() if g is not None])
        chrom = match.groups()[0]
        if num_groups == 1:
            return 'not_found', {'message': 'The search for {} returned no results.'.format(query)}
        pos = int(match.groups()[1].replace(',',''))
        if num_groups == 2:
            return 'region', {'chrom': chrom, 'start':pos, 'stop':pos}
        if num_groups == 3:
            end = int(match.groups()[2].replace(',',''))
            return 'region', {'chrom': chrom, 'start':pos, 'stop':end}
        return 'variant', {'variant_id': '{}-{}-{}-{}'.format(chrom, pos, match.groups()[2], match.groups()[3])}
    # return 'not_found', {'query': query_orig}
    return 'not_found', {'message': 'The search for {} returned no results.'.format(query)}


class IntervalSet(object):
    EXON_PADDING = 20

    def __init__(self, chrom, list_of_pairs):
        self.chrom = chrom
        self._list_of_pairs = list_of_pairs # [[start1, stop1], [start2, stop2], ...]
    @classmethod
    def from_chrom_start_stop(cls, chrom, start, stop):
        Xpos.check_chrom(chrom)
        assert start <= stop
        return cls(chrom, [[start, stop]])
    @classmethod
    def from_xstart_xstop(cls, xstart, xstop):
        chrom1, start = Xpos.to_chrom_pos(xstart)
        chrom2, stop = Xpos.to_chrom_pos(xstop)
        assert start <= stop
        assert chrom1 == chrom2
        return cls(chrom1, [[start, stop]])
    @classmethod
    def from_gene(cls, db, gene_id):
        exons = db.exons.find({'gene_id': gene_id, 'feature_type': { "$in": ['CDS', 'UTR', 'exon'] }}, projection={'_id': False})
        return cls._from_exons(exons)
    @classmethod
    def from_transcript(cls, db, transcript_id):
        exons = db.exons.find({'transcript_id': transcript_id, 'feature_type': { "$in": ['CDS', 'UTR', 'exon'] }}, projection={'_id': False})
        return cls._from_exons(exons)
    @classmethod
    def _from_exons(cls, exons):
        # note: these "exons" are not all literally exons, some are CDS or UTR features
        exons = sorted(list(exons), key=lambda exon: exon['start'])
        assert len(exons) > 0
        try:
            assert boltons.iterutils.same(exon['chrom'] for exon in exons)
        except AssertionError:
            # if we're in PAR, only show X.  TODO: actually solve this for the long term.
            assert all(exon['chrom'] in ['X','Y'] for exon in exons)
            exons = [exon for exon in exons if exon['chrom'] == 'X']
        regions = []
        for exon in exons:
            assert exon['start'] <= exon['stop'] # There are some exons with start==stop, which I don't understand
            start, stop = exon['start']-cls.EXON_PADDING, exon['stop']+cls.EXON_PADDING
            if not regions or regions[-1][1] <= start:
                regions.append([start, stop])
            elif regions[-1][1] < stop:
                regions[-1][1] = stop
        return cls(exons[0]['chrom'], regions)

    def to_obj(self):
        return {'chrom': self.chrom, 'list_of_pairs': self._list_of_pairs}
    def to_mongo(self):
        return {'$or': self.to_list_of_mongos()}
    def to_list_of_mongos(self):
        return [{'xpos': {'$gte':Xpos.from_chrom_pos(self.chrom,start),'$lte':Xpos.from_chrom_pos(self.chrom,stop)}} for (start,stop) in self._list_of_pairs]
    def __str__(self):
        return '{}:{}'.format(self.chrom, ','.join('{}-{}'.format(*pair) for pair in self._list_of_pairs))
    __repr__ = __str__

    def get_start(self): return self._list_of_pairs[0][0]
    def get_stop(self): return self._list_of_pairs[-1][1]
    def get_length(self): return sum(pair[1] - pair[0] for pair in self._list_of_pairs)
    def to_region_dict(self): return {'chrom': self.chrom, 'start':self.get_start(), 'stop':self.get_stop()}
    def to_region_dashed(self): return '{}-{}-{}'.format(self.chrom, self.get_start(), self.get_stop())


class TranscriptSet(object):
    # TODO: maybe just make each of these methods return the json?  or will this class be more complex?
    def __init__(self, genes):
        self.genes = genes
    @classmethod
    def from_gene(cls, db, gene_id):
        all_exons = list(db.exons.find({'gene_id': gene_id}, {'_id':False}))
        return cls._from_exons(db, all_exons)
    @classmethod
    def from_transcript(cls, db, transcript_id):
        all_exons = list(db.exons.find({'transcript_id': transcript_id}, {'_id':False}))
        return cls._from_exons(db, all_exons)
    @classmethod
    def from_chrom_start_stop(cls, db, chrom, start, stop):
        xstart,xstop = Xpos.from_chrom_pos(chrom,start),Xpos.from_chrom_pos(chrom,stop)
        all_exons = list(db.exons.find({'xstop':{'$gte':xstart},'xstart':{'$lte':xstop}}, {'_id':False}))
        return cls._from_exons(db, all_exons)
    @classmethod
    def _from_exons(cls, db, all_exons):
        '''return is like [{gene_name:'PCSK9', gene_id:'ENSG123', transcripts:[{transcript_id:'ENST234',start,stop,exons:[{start,stop,strand,feature_type}]}]}]'''
        for exon in all_exons: assert exon['feature_type'] in ['exon', 'CDS', 'UTR'] and exon['strand'] in ['+','-']
        all_transcripts = []
        for transcript_id, exons in sortedgroupby(all_exons, key=lambda exon:exon['transcript_id']):
            exons = sorted(exons, key=lambda exon:exon['start'])
            transcript = db.transcripts.find_one({'transcript_id':transcript_id})
            gene_id = exons[0]['gene_id']
            exons = [{key: exon[key] for key in ['feature_type','strand','start','stop']} for exon in exons]
            weight = 0 # 10 * CDS length + UTR length + exon length + 1e10 * canonical
            for exon in exons:
                length = exon['stop']+1 - exon['start']
                weight += length * {'CDS':10, 'UTR':1, 'exon':1}[exon['feature_type']]
            all_transcripts.append({
                'gene_id':gene_id,'transcript_id':transcript_id,
                'start':transcript['start'],'stop':transcript['stop'],
                'exons':exons,'weight':weight
            })
        genes = []
        for gene_id, transcripts in sortedgroupby(all_transcripts, key=lambda trans:trans['gene_id']):
            gene = get_gene(db, gene_id)
            gene_name = gene['gene_name'] if gene else None
            canonical_transcript_id = gene.get('canonical_transcript') if gene else None
            transcripts = list(transcripts)
            for transcript in transcripts:
                if transcript['transcript_id'] == canonical_transcript_id:
                    transcript['canonical'] = True
                    transcript['weight'] += 1e10
            transcripts = sorted(transcripts, key=lambda trans:-trans['weight'])
            genes.append({
                'gene_id': gene_id,
                'gene_name': gene_name,
                'transcripts': transcripts,
            })
        genes.sort(key=lambda gene:-gene['transcripts'][0]['weight'])
        for gene in genes:
            for transcript in gene['transcripts']:
                del transcript['gene_id']
                del transcript['weight']
        return cls(genes)


def get_metrics(db):
    metrics = []
    cursor = db.metrics.find({'type': 'percentiles'}, projection = {'_id': False})
    for document in cursor:
        metrics.append(document)
    return metrics


def remove_some_extraneous_information(variant):
    """Remove information not needed by variant.html or any other page"""
    for key in ['xpos','xstop','vep_annotations',]: variant.pop(key, None)



def get_summary_for_intervalset(db, intervalset):
    # Note: querying for each extent in intervalset.to_list_of_mongos() is >100X faster than using intervalset.to_mongo() and I have no idea why. Try query planner?
    st = time.time()
    mongo_match_cond = {
        'lof': {'$lt': ['$worst_csqidx', Consequence.as_obj['n_lof']]},
        'lof_lc': {'$and': [
            {'$lt': ['$worst_csqidx', Consequence.as_obj['n_lof']]},
            {'$in':[{'k':"LoF",'v':"LC"},{'$objectToArray':{'$arrayElemAt':["$vep_annotations",0]}}]} # this is gross but I couldn't find a better way.
        ]},
        'mis': {'$and': [{'$gte': ['$worst_csqidx', Consequence.as_obj['n_lof']]},     {'$lt':['$worst_csqidx', Consequence.as_obj['n_lof_mis']]}]},
        'syn': {'$and': [{'$gte': ['$worst_csqidx', Consequence.as_obj['n_lof_mis']]}, {'$lt':['$worst_csqidx', Consequence.as_obj['n_lof_mis_syn']]}]},
        'indel': {'$or': [{'$ne': [1, {'$strLenBytes':'$ref'}]}, {'$ne': [1, {'$strLenBytes':'$alt'}]}]},
    }
    keys = 'lof lof_lc mis syn indel total'.split()
    ret = {key:0 for key in keys}
    for mongo_match_region in intervalset.to_list_of_mongos():
        x = db.variants.aggregate([
            {'$match': mkdict(mongo_match_region, {'filter':'PASS'})},
            {'$group': {
                '_id': None,
                'lof':  {'$sum':{'$cond':[mongo_match_cond['lof'],1,0]}},
                'lof_lc':  {'$sum':{'$cond':[mongo_match_cond['lof_lc'],1,0]}},
                'mis':  {'$sum':{'$cond':[mongo_match_cond['mis'],1,0]}},
                'syn':  {'$sum':{'$cond':[mongo_match_cond['syn'],1,0]}},
                'indel':{'$sum':{'$cond':[mongo_match_cond['indel'],1,0]}},
                'total':{'$sum':1},
            }},
        ])
        x = list(x);
        if len(x) == 0: continue # no variants in interval
        assert len(x) == 1; x = x[0]
        for key in keys: ret[key] += x.get(key,0)
    print '## SUMMARY: spent {:0.3f} seconds tabulating {} variants'.format(time.time() - st, ret['total'])
    return [
        ('All - SNPs', ret['total'] - ret['indel']),
        ('All - Indels', ret['indel']),
        ('Coding - LoF', ret['lof']),
        ('Coding - LoF - Low Confidence', ret['lof_lc']),
        ('Coding - Missense', ret['mis']),
        ('Coding - Synonymous', ret['syn']),
    ]


def get_variants_subset_for_intervalset(db, intervalset, columns_to_return, order, filter_info, skip, length):
    # 1. match what the user asked for - using [intervalset, filter_info]
    # 2. project to just keys for sorting, sort, and get `_id`s - using [order]
    # 3. get `n_filtered` and `length`-many `_id`s - using [skip, length]
    # 4. look up those `_id`s and project - using [columns_to_return]
    st = time.time()

    mongo_match = [intervalset.to_mongo()]
    if filter_info.get('filter_value',None) is not None:
        if filter_info['filter_value'] == 'PASS': mongo_match.append({'filter': 'PASS'})
        elif filter_info['filter_value'] == 'not PASS': mongo_match.append({'filter': {'$ne': 'PASS'}})
    if isinstance(filter_info.get('maf_ge',None),(float,int)):
        assert 0 <= filter_info['maf_ge'] <= 0.5
        if filter_info['maf_ge'] > 0: mongo_match.append({'$and': [{'allele_freq': {'$gte': filter_info['maf_ge']}},{'allele_freq': {'$lte': 1-filter_info['maf_ge']}}]})
    if isinstance(filter_info.get('maf_le',None),(float,int)):
        assert 0 <= filter_info['maf_le'] <= 0.5
        if filter_info['maf_le'] < 0.5: mongo_match.append({'$or': [{'allele_freq': {'$lte': filter_info['maf_le']}},{'allele_freq': {'$gte': 1-filter_info['maf_le']}}]})
    if filter_info.get('category',None) is not None:
        if filter_info['category'].strip() == 'LoF': mongo_match.append({'worst_csqidx': {'$lt': Consequence.as_obj['n_lof']}})
        elif filter_info['category'].strip() == 'LoF+Missense': mongo_match.append({'worst_csqidx': {'$lt': Consequence.as_obj['n_lof_mis']}})

    cols = {
        # after pre-processing, these will look like:
        # <name>: {'sort': {'project': <projection>, 'sort_key': <key>}, 'return': {'project': <projection>}}
        # <name>: {'sort': False, 'return': {'project': <projection>}}
        'allele': {'return': ['rsids', 'ref', 'alt']},
        'pos': {'sort': 'xpos'},
        'csq': {'sort': 'worst_csqidx', 'return':{'project': {
            'worst_csqidx':1,
            'HGVS':'$worst_csq_HGVS',
            'low_conf': {'$in':[{'k':"LoF",'v':"LC"},{'$objectToArray':{'$arrayElemAt':["$vep_annotations",0]}}]}, # Gross, but I don't know a better way.
        }}},
        'filter': {},
        'allele_count': {'sort': True},
        'allele_num': {'sort': True},
        'het': {'sort': {'project': {'het': {'$subtract':['$allele_count',{'$multiply':[2,'$hom_count']}]}}, 'sort_key': 'het'},
                'return': {'project': {'het': {'$subtract':['$allele_count',{'$multiply':[2,'$hom_count']}]}}}},
        'hom_count': {'sort': True},
        'allele_freq': {'sort': True},
        'cadd_phred': {'sort': True},
    }
    for name, col in cols.items():
        try:
            if 'sort' not in col: col['sort'] = False
            if col['sort'] == True: col['sort'] = name
            if isinstance(col['sort'], str): col['sort'] = {'project': {col['sort']:1}, 'sort_key':col['sort']}
            assert col['sort'] == False or isinstance(col['sort']['project'], dict) and isinstance(col['sort']['sort_key'], str)
            if 'return' not in col: col['return'] = [name]
            if isinstance(col['return'], list): col['return'] = {'project': {k:1 for k in col['return']}}
            assert isinstance(col['return']['project'], dict)
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

    mongo_projection = mkdict(*[cols[ctr['name']]['return']['project'] for ctr in columns_to_return], _id=False)

    v_ids_curs = db.variants.aggregate([
        {'$match': {'$and': mongo_match}},
        {'$project': mongo_projection_before_sort},
        {'$sort': mongo_sort},
        {'$project': {'_id': 1}},
        {'$group': {'_id':0, 'count':{'$sum':1}, 'results':{'$push':'$$ROOT'}}},
        {'$project': {'_id':0, 'count':1, 'ids':{'$slice':['$results',skip,length]}}},
    ])
    print '## VARIANT_SUBSET: spent {:.3f} seconds creating cursor'.format(time.time()-st); st = time.time()
    v_ids_result = list(v_ids_curs)
    if len(v_ids_result) == 0:
        n_filtered, variants = 0, []
    else:
        assert len(v_ids_result) == 1
        n_filtered = v_ids_result[0]['count']
        print '## VARIANT_SUBSET: spent {:0.3f} seconds counting {} variants that match filters'.format(time.time()-st, n_filtered); st = time.time()
        v_ids = [v['_id'] for v in v_ids_result[0]['ids']]
        variants = [next(db.variants.aggregate([{'$match': {'_id': vid}}, {'$project': mongo_projection}])) for vid in v_ids] # b/c fancy projections require .aggregate()
        print '## VARIANT_SUBSET: spent {:0.3f} seconds fetching {} full variants by id'.format(time.time()-st, len(variants)); st = time.time()

    return {
        'recordsFiltered': n_filtered,
        'recordsTotal': n_filtered,
        'data': variants
    }


def get_variants_csv_str_for_intervalset(db, intervalset):
    import io, csv
    out = io.BytesIO()
    writer = csv.writer(out)
    fields = 'chrom pos ref alt rsids filter genes allele_num allele_count allele_freq hom_count site_quality quality_metrics.DP cadd_phred'.split()
    writer.writerow(fields)
    variants = get_variants_in_intervalset(db, intervalset)
    for v in variants:
        row = []
        for field in fields:
            if '.' in field: parts = field.split('.', 1); row.append(v.get(parts[0], {}).get(parts[1], ''))
            elif field in ['rsids','genes']: row.append('|'.join(v.get(field, [])))
            else: row.append(v.get(field, ''))
        writer.writerow(row)
    return out.getvalue()
def get_variants_in_intervalset(db, intervalset):
    """Variants that overlap an intervalset"""
    for mongo_match_region in intervalset.to_list_of_mongos():
        for variant in db.variants.find(mongo_match_region, projection={'_id': False}):
            yield variant
