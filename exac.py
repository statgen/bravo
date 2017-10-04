#!/usr/bin/env python2

import itertools
import json
import os
import pymongo
import pysam
import gzip
import random
import boltons.cacheutils

from flask import Flask, Response, request, session, g, redirect, url_for, abort, render_template, flash, jsonify, make_response, send_file, Blueprint
from flask_compress import Compress
from flask_errormail import mail_on_500
from flask_login import LoginManager, UserMixin, login_user, logout_user, current_user

from collections import defaultdict, Counter
from multiprocessing import Process
import multiprocessing
import glob
import traceback
import time
import sys
import functools
import contextlib

from parsing import *
import lookups
from lookups import IntervalSet, TranscriptSet
from utils import *
from base_coverage import CoverageHandler
import auth

bp = Blueprint('bp', __name__, template_folder='templates', static_folder='static')

app = Flask(__name__)
app.config.from_object('flask_config.BravoFreeze5GRCh38Config')
mail_on_500(app, app.config['ADMINS'])
app.config['COMPRESS_LEVEL'] = 2 # Since we don't cache, faster=better
Compress(app)
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 5 # 5 second browser cache timeout
app.config['TEMPLATES_AUTO_RELOAD'] = True

MAX_REGION_LENGTH = int(350e3) # Longer than TTN (305kb), short enough to perform okay.

def get_db(new_connection=False):
    # Only use the database within a request context! Something about threads/forks.
    # See <https://jira.mongodb.org/browse/PYTHON-961>
    # Note: I just added `connect=False`, so maybe we don't need this function anymore (unless used with new_connection=True)
    if new_connection:
        client = pymongo.MongoClient(host=app.config['MONGO']['host'], port=app.config['MONGO']['port'], connect=False)
    else:
        client = get_db._mongo_client
    return client[app.config['MONGO']['name']]
get_db._mongo_client = pymongo.MongoClient(host=app.config['MONGO']['host'], port=app.config['MONGO']['port'], connect=False)

@boltons.cacheutils.cached({})
def get_autocomplete_strings():
    autocomplete_strings = get_db().genes.distinct('gene_name')
    autocomplete_strings.extend(get_db().genes.distinct('other_names', {'other_names': {'$ne': None}}))
    return sorted(set(autocomplete_strings))

@boltons.cacheutils.cached({})
def get_coverage_handler():
    return CoverageHandler(app.config['BASE_COVERAGE'])


def get_tabix_file_contig_pairs(tabix_filenames):
    filename_contig_pairs = []
    for tabix_filename in tabix_filenames:
        with pysam.Tabixfile(tabix_filename) as tabix_file:
            for contig in tabix_file.contigs:
                filename_contig_pairs.append((tabix_filename, contig))
    def _sort_key(pair):
        (filename, contig) = pair
        if contig.startswith('chr'): contig = contig[3:]
        if contig.isdigit(): return int(contig)
        return 0 # for X/Y/MT
    filename_contig_pairs.sort(key=_sort_key) # Sort from large -> small chromosomes
    return filename_contig_pairs

def get_records_from_tabix_contig(tabix_filename, contig, record_parser):
    start_time = time.time()
    with pysam.Tabixfile(tabix_filename) as tabix_file:
        record_i = 0 # in case record_parser never yields anything.
        for record_i, parsed_record in enumerate(record_parser(itertools.chain(tabix_file.header, tabix_file.fetch(contig, 0, 10**10, multiple_iterators=True))), start=1):
            yield parsed_record

            if record_i % int(1e6) == 0:
                print("Loaded {:11,} records in {:6,} seconds from contig {!r:6} of {!r}".format(record_i, int(time.time()-start_time), contig, tabix_filename))
    print("Loaded {:11,} records in {:6,} seconds from contig {!r:6} of {!r}".format(record_i, int(time.time()-start_time), contig, tabix_filename))

def _load_variants_from_tabix_file_and_contig(args):
    tabix_file, contig = args
    db = get_db(new_connection=True)
    variants_generator = get_records_from_tabix_contig(tabix_file, contig, get_variants_from_sites_vcf)
    try:
        db.variants.insert(variants_generator, w=0)
    except pymongo.errors.InvalidOperation:
        pass  # handle error when variant_generator is empty

def load_variants_file():
    db = get_db()
    db.variants.drop()
    print("Dropped db.variants")

    if len(app.config['SITES_VCFS']) == 0:
        raise IOError("No vcf file found")

    with contextlib.closing(multiprocessing.Pool(app.config['LOAD_DB_PARALLEL_PROCESSES'])) as pool:
        # workaround for Pool.map() from <http://stackoverflow.com/a/1408476/1166306>
        pool.map_async(_load_variants_from_tabix_file_and_contig, get_tabix_file_contig_pairs(app.config['SITES_VCFS'])).get(9999999)

    # TODO: use db.variants.create_indexes([pymongo.operations.IndexModel(key) for key in 'xpos xstop rsids filter'.split()])
    for key in 'xpos xstop rsids filter'.split():
        print 'creating index on', key
        db.variants.create_index(key)


def load_gene_models():
    db = get_db()

    db.genes.drop()
    db.transcripts.drop()
    db.exons.drop()
    print 'Dropped db.genes, db.transcripts, and db.exons.'

    start_time = time.time()

    canonical_transcripts = {}
    with gzip.open(app.config['CANONICAL_TRANSCRIPT_FILE']) as canonical_transcript_file:
        for gene, transcript in get_canonical_transcripts(canonical_transcript_file):
            canonical_transcripts[gene] = transcript

    omim_annotations = {}
    with gzip.open(app.config['OMIM_FILE']) as omim_file:
        for fields in get_omim_associations(omim_file):
            if fields is None:
                continue
            gene, transcript, accession, description = fields
            omim_annotations[gene] = (accession, description)

    dbnsfp_info = {}
    with gzip.open(app.config['DBNSFP_FILE']) as dbnsfp_file:
        for dbnsfp_gene in get_dbnsfp_info(dbnsfp_file):
            other_names = [other_name.upper() for other_name in dbnsfp_gene['gene_other_names']]
            dbnsfp_info[dbnsfp_gene['ensembl_gene']] = (dbnsfp_gene['gene_full_name'], other_names)

    print 'Done loading metadata. Took %s seconds' % int(time.time() - start_time)

    # grab genes from GTF
    start_time = time.time()
    with gzip.open(app.config['GENCODE_GTF']) as gtf_file:
        for gene in get_genes_from_gencode_gtf(gtf_file):
            gene_id = gene['gene_id']
            if gene_id in canonical_transcripts:
                gene['canonical_transcript'] = canonical_transcripts[gene_id]
            if gene_id in omim_annotations:
                gene['omim_accession'] = omim_annotations[gene_id][0]
                gene['omim_description'] = omim_annotations[gene_id][1]
            if gene_id in dbnsfp_info:
                gene['full_gene_name'] = dbnsfp_info[gene_id][0]
                gene['other_names'] = dbnsfp_info[gene_id][1]
            db.genes.insert(gene, w=0)

    print 'Done loading genes. Took %s seconds' % int(time.time() - start_time)

    start_time = time.time()
    db.genes.ensure_index('gene_id')
    db.genes.ensure_index('gene_name_upper')
    db.genes.ensure_index('gene_name')
    db.genes.ensure_index('other_names')
    db.genes.ensure_index('xstart')
    db.genes.ensure_index('xstop')
    print 'Done indexing gene table. Took %s seconds' % int(time.time() - start_time)

    # and now transcripts
    start_time = time.time()
    with gzip.open(app.config['GENCODE_GTF']) as gtf_file:
        db.transcripts.insert((transcript for transcript in get_transcripts_from_gencode_gtf(gtf_file)), w=0)
    print 'Done loading transcripts. Took %s seconds' % int(time.time() - start_time)

    start_time = time.time()
    db.transcripts.ensure_index('transcript_id')
    db.transcripts.ensure_index('gene_id')
    print 'Done indexing transcript table. Took %s seconds' % int(time.time() - start_time)

    # Building up gene definitions
    start_time = time.time()
    with gzip.open(app.config['GENCODE_GTF']) as gtf_file:
        db.exons.insert((exon for exon in get_exons_from_gencode_gtf(gtf_file)), w=0)
    print 'Done loading exons. Took %s seconds' % int(time.time() - start_time)

    start_time = time.time()
    db.exons.ensure_index('exon_id')
    db.exons.ensure_index('transcript_id')
    db.exons.ensure_index('gene_id')
    print 'Done indexing exon table. Took %s seconds' % int(time.time() - start_time)


def _load_dbsnp_from_tabix_file_and_contig(args):
    dbsnp_file, contig = args
    db = get_db(new_connection=True)
    dbsnp_record_generator = get_records_from_tabix_contig(dbsnp_file, contig, get_snp_from_dbsnp_file)
    try:
        db.dbsnp.insert(dbsnp_record_generator, w=0)
    except pymongo.errors.InvalidOperation:
        pass  # handle error when generator is empty

def load_dbsnp_file():
    db = get_db()

    db.dbsnp.drop()
    db.dbsnp.ensure_index('rsid') # It seems faster to build these indexes before inserts.  Strange.
    db.dbsnp.ensure_index('xpos')
    start_time = time.time()
    dbsnp_file = app.config['DBSNP_FILE']

    print "Loading dbsnp from %s" % dbsnp_file
    if os.path.isfile(dbsnp_file + ".tbi"):
        with contextlib.closing(multiprocessing.Pool(app.config['LOAD_DB_PARALLEL_PROCESSES'])) as pool:
            # workaround for Pool.map() from <http://stackoverflow.com/a/1408476/1166306>
            pool.map_async(_load_dbsnp_from_tabix_file_and_contig, get_tabix_file_contig_pairs([dbsnp_file])).get(9999999)
        print('Done loading dbSNP in {:,} seconds'.format(int(time.time() - start_time)))

    elif os.path.isfile(dbsnp_file):
        # see if non-tabixed .gz version exists
        print(("WARNING: %(dbsnp_file)s.tbi index file not found. Will use single thread to load dbsnp."
               "To create a tabix-indexed dbsnp file based on UCSC dbsnp, do: \n"
               "   wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp141.txt.gz \n"
               "   gzcat snp141.txt.gz | cut -f 1-5 | bgzip -c > snp141.txt.bgz \n"
               "   tabix -0 -s 2 -b 3 -e 4 snp141.txt.bgz") % locals())
        with gzip.open(dbsnp_file) as f:
            db.dbsnp.insert((snp for snp in get_snp_from_dbsnp_file(f)), w=0)

    else:
        raise Exception("dbsnp file %s(dbsnp_file)s not found." % locals())

def precalculate_metrics():
    import numpy
    db = get_db()
    total_num_variants = db.variants.count()
    # Mongo sorts the whole table if I request >5% of all documents, so I take only 4% for safety.
    num_variants_to_sample = clamp(total_num_variants/25, min_value=int(10e3), max_value=min(int(1e6),total_num_variants))
    print 'Sampling {:,} of {:,} variants...'.format(num_variants_to_sample, total_num_variants)
    metrics = defaultdict(Counter)
    qualities_by_af = defaultdict(Counter)
    start_time = time.time()
    variant_iterator = db.variants.aggregate([
        {'$sample': {'size':num_variants_to_sample}},
        {'$project': {'quality_metrics':1, 'site_quality':1, 'allele_num':1, 'allele_count':1, '_id':0}},
    ])
    for variant_i, variant in enumerate(variant_iterator):
        if variant_i % int(100e3) == 0:
            num_keys = sum(len(cntr.viewkeys()) for cntr in metrics.values())
            print 'Read {:,} variants. Took {:.0f} seconds. Made {:,} distinct keys'.format(variant_i, time.time() - start_time, num_keys)
            start_time = time.time()
        if 'DP' in variant['quality_metrics'] and float(variant['quality_metrics']['DP']) == 0:
            print('Warning: variant with id {} has depth of 0'.format(variant['_id']))
        for metric, value in variant['quality_metrics'].iteritems():
            metrics[metric][float(value)] += 1
        qual = float(variant['site_quality'])
        metrics['site_quality'][qual] += 1
        if variant['allele_num'] == 0: continue
        if variant['allele_count'] == 1:
            qualities_by_af['singleton'][qual] += 1
        elif variant['allele_count'] == 2:
            qualities_by_af['doubleton'][qual] += 1
        else:
            variant_af = float(variant['allele_count'])/variant['allele_num']
            for bucket_af in AF_BUCKETS:
                if variant_af < bucket_af:
                    qualities_by_af[bucket_af][qual] += 1
                    break
    num_keys = sum(len(cntr.viewkeys()) for cntr in metrics.values())
    print 'Read {:,} variants. Took {:.0f} seconds. Made {:,} distinct keys'.format(variant_i, time.time() - start_time, num_keys)
    print 'Done reading variants. Dropping metrics database... '
    db.metrics.drop()
    start_time = time.time()
    print 'Dropped metrics database. Calculating metrics...'
    for metric in metrics:
        print '... for metric', metric
        bin_range = None
        data = metrics[metric]
        if metric == 'DP': data = {numpy.log(key): count for key,count in data.iteritems()}
        if metric == 'FS': bin_range = (0, 20)
        elif metric == 'VQSLOD': bin_range = (-20, 20)
        elif metric == 'InbreedingCoeff': bin_range = (0, 1)
        hist = histogram_from_counter(data, num_bins=40, bin_range=bin_range)
        db.metrics.insert({'metric': metric, 'mids': hist['left_edges'], 'hist': hist['counts']})
    for af_bin in qualities_by_af:
        print '... for site quality for af', af_bin
        data = {numpy.log(key): count for key,count in qualities_by_af[af_bin].iteritems()}
        hist = histogram_from_counter(data, num_bins=40)
        db.metrics.insert({'metric': 'binned_%s' % af_bin, 'mids': hist['mids'], 'hist': hist['counts']})
    db.metrics.ensure_index('metric')
    print '... (took {:.0f} seconds)'.format(time.time() - start_time)
    print 'Done pre-calculating metrics!'

def create_users():
    db = get_db()
    db.users.drop()
    print 'Dropped users database.'
    db.users.ensure_index('user_id')
    print 'Created new users database.'


def require_agreement_to_terms_and_store_destination(func):
    """
    This decorator for routes checks that the user is logged in and has agreed to the terms.
    If they haven't, their intended destination is stored and they're sent to get authorized.
    I think that it has to be placed AFTER @app.route() so that it can capture `request.path`.
    """
    # inspired by <https://flask-login.readthedocs.org/en/latest/_modules/flask_login.html#login_required>
    @functools.wraps(func)
    def decorated_view(*args, **kwargs):
        if hasattr(current_user, 'agreed_to_terms') and current_user.agreed_to_terms:
            return func(*args, **kwargs)
        else:
            print('unauthorized user {!r} visited the url [{!r}]'.format(current_user, request.path))
            session['original_destination'] = request.path
            return redirect(url_for('.get_authorized'))
        return func(*args, **kwargs)
    return decorated_view

def _log(message=''):
    url = request.full_path.rstrip('?')
    if url.startswith(app.config['URL_PREFIX']): url = url[len(app.config['URL_PREFIX']):]
    print('{}  {}{}'.format(current_user, url, message))

def _err():
    url = request.full_path.rstrip('?')
    if url.startswith(app.config['URL_PREFIX']): url = url[len(app.config['URL_PREFIX']):]
    error = traceback.format_exc()
    if request.form: print('Failed on {} with form {} and error:\n{}'.format(url, request.form, error))
    else: print('Failed on {} with error:\n{}'.format(url, error))

@bp.route('/')
def homepage():
    return render_template('homepage.html')


@bp.route('/api/autocomplete')
def autocomplete():
    query = request.args.get('query', '')
    suggestions = lookups.get_awesomebar_suggestions(get_autocomplete_strings(), query)
    _log('  =>  {} results'.format(len(suggestions)))
    return jsonify([{'value': s} for s in sorted(suggestions)])


@bp.route('/awesome')
def awesome():
    db = get_db()
    query = request.args.get('query')
    datatype, redirect_args = lookups.get_awesomebar_result(db, query)
    _log('  =>  {}_page({})'.format(datatype, redirect_args))
    return redirect(url_for('.{}_page'.format(datatype), **redirect_args))


@bp.route('/profile', methods=['GET', 'POST'])
@require_agreement_to_terms_and_store_destination
def user_profile_page():
    try:
        _log()
        error = None
        success = None
        if request.method == 'POST':
            enabled_api = False if request.form.get('enabled_api', '').lower() != 'on' else True
            no_newsletters = False if request.form.get('no_newsletters', '').lower() != 'on' else True
            if error is None:
                db = get_db()
                result = db.users.update_one({"user_id": current_user.get_id()}, {"$set": {"enabled_api": enabled_api, "no_newsletters": no_newsletters}})
                success = True
                current_user.enabled_api = enabled_api
                current_user.no_newsletters = no_newsletters
        return render_template('user_profile.html', error = error, success = success)
    except: _err(); abort(404)

@bp.route('/variant/<variant_id>')
@require_agreement_to_terms_and_store_destination
def variant_page(variant_id):
    db = get_db()
    try:
        _log()
        variant = lookups.get_variant_by_variant_id(db, variant_id, default_to_boring_variant=True)
        if not variant: return error_page('Variant {!r} not found'.format(variant_id))
        if 'pop_afs' in variant: variant['pop_afs'][app.config['DATASET_NAME']] = variant['allele_freq']

        consequence_drilldown = ConsequenceDrilldown.from_variant(variant)
        gene_for_top_csq, top_HGVSs = ConsequenceDrilldown.get_top_gene_and_HGVSs(consequence_drilldown)
        consequence_drilldown_columns = ConsequenceDrilldown.split_into_two_columns(consequence_drilldown)

        base_coverage = get_coverage_handler().get_coverage_for_intervalset(
            IntervalSet.from_xstart_xstop(variant['xpos'], variant['xpos']+len(variant['ref'])-1))
        metrics = lookups.get_metrics(db, variant)
        lookups.remove_some_extraneous_information(variant)

        return render_template(
            'variant.html',
            variant=variant,
            base_coverage=base_coverage,
            consequences=consequence_drilldown,
            consequence_columns=consequence_drilldown_columns,
            any_covered=bool(base_coverage),
            metrics=metrics,
            top_HGVSs=top_HGVSs,
            gene_for_top_csq=gene_for_top_csq,
        )
    except: _err(); abort(404)


@bp.route('/gene/<gene_id>')
@require_agreement_to_terms_and_store_destination
def gene_page(gene_id):
    db = get_db()
    try:
        gene = lookups.get_gene(db, gene_id)
        if not gene: return error_page('Gene {!r} not found'.format(gene_id))
        _log('   ({})'.format(gene.get('gene_name')))
        intervalset = IntervalSet.from_gene(db, gene_id)
        genes = TranscriptSet.from_gene(db, gene_id).genes
        return render_template(
            'gene.html',
            intervalset=intervalset, genes=genes, csq=Consequence.as_obj,
            gene=gene,
        )
    except:_err(); abort(404)

@bp.route('/transcript/<transcript_id>')
@require_agreement_to_terms_and_store_destination
def transcript_page(transcript_id):
    db = get_db()
    try:
        _log()
        transcript = lookups.get_transcript(db, transcript_id)
        if not transcript: return error_page('Transcript {!r} not found'.format(transcript_id))
        gene = lookups.get_gene(db, transcript['gene_id'])
        intervalset = IntervalSet.from_transcript(db, transcript_id)
        genes = TranscriptSet.from_transcript(db, transcript_id).genes
        return render_template(
            'transcript.html',
            intervalset=intervalset, genes=genes, csq=Consequence.as_obj,
            gene=gene,
            transcript=transcript,
        )
    except:_err(); abort(404)

@bp.route('/region/<chrom>-<start>-<stop>')
@require_agreement_to_terms_and_store_destination
def region_page(chrom, start, stop):
    db = get_db()
    try:
        _log()
        try: start,stop = int(start),int(stop)
        except: return error_page('Positions not integers: {!r}, {!r}'.format(start, stop))
        if start > stop: return error_page("The region '{chrom}-{start}-{stop}' stops before it starts. Did you mean '{chrom}-{stop}-{start}'?".format(chrom=chrom, start=start, stop=stop))
        if stop-start > MAX_REGION_LENGTH: return error_page("The region '{chrom}-{start}-{stop}' is {:,} bases. We only accept regions shorter than {:,} bases.".format(stop-start, MAX_REGION_LENGTH, chrom=chrom, start=start, stop=stop))
        if start == stop: start -= 20; stop += 20

        intervalset = IntervalSet.from_chrom_start_stop(chrom, start, stop)
        genes = TranscriptSet.from_chrom_start_stop(db, chrom, start, stop).genes
        return render_template(
            'region.html',
            intervalset=intervalset, genes=genes, csq=Consequence.as_obj,
        )
    except: _err(); abort(404)


@bp.route('/download/gene/<gene_id>')
@require_agreement_to_terms_and_store_destination
def download_gene_variants(gene_id):
    db = get_db()
    try:
        intervalset = IntervalSet.from_gene(get_db(), gene_id)
        return _get_variants_csv_for_intervalset(intervalset, '{}.csv'.format(gene_id))
    except:_err(); abort(404)
@bp.route('/download/transcript/<transcript_id>')
@require_agreement_to_terms_and_store_destination
def download_transcript_variants(transcript_id):
    db = get_db()
    try:
        intervalset = IntervalSet.from_transcript(get_db(), transcript_id)
        return _get_variants_csv_for_intervalset(intervalset, '{}.csv'.format(transcript_id))
    except:_err(); abort(404)
@bp.route('/download/region/<chrom>-<start>-<stop>')
@require_agreement_to_terms_and_store_destination
def download_region_variants(chrom, start, stop):
    try:
        start,stop = int(start),int(stop); assert stop-start <= MAX_REGION_LENGTH
        intervalset = IntervalSet.from_chrom_start_stop(chrom, start, stop)
        return _get_variants_csv_for_intervalset(intervalset, 'chr{}-{}-{}.csv'.format(chrom, start, stop))
    except:_err(); abort(404)
def _get_variants_csv_for_intervalset(intervalset, filename):
    _log()
    resp = make_response(lookups.get_variants_csv_str_for_intervalset(get_db(), intervalset))
    resp.headers['Content-Disposition'] = 'attachment; filename={}'.format(filename)
    resp.mimetype='text/csv'
    return resp


@bp.route('/api/summary/gene/<gene_id>')
@require_agreement_to_terms_and_store_destination
def gene_summary_api(gene_id):
    try:
        intervalset = IntervalSet.from_gene(get_db(), gene_id)
        return jsonify(lookups.get_summary_for_intervalset(get_db(), intervalset))
    except:_err(); abort(404)
@bp.route('/api/summary/transcript/<transcript_id>')
@require_agreement_to_terms_and_store_destination
def transcript_summary_api(transcript_id):
    try:
        intervalset = IntervalSet.from_transcript(get_db(), transcript_id)
        return jsonify(lookups.get_summary_for_intervalset(get_db(), intervalset))
    except:_err(); abort(404)
@bp.route('/api/summary/region/<chrom>-<start>-<stop>')
@require_agreement_to_terms_and_store_destination
def region_summary_api(chrom, start, stop):
    try:
        start,stop = int(start),int(stop); assert stop-start <= MAX_REGION_LENGTH
        intervalset = IntervalSet.from_chrom_start_stop(chrom, start, stop)
        return jsonify(lookups.get_summary_for_intervalset(get_db(), intervalset))
    except:_err(); abort(404)


@bp.route('/api/variants/gene/<gene_id>', methods=['POST'])
@require_agreement_to_terms_and_store_destination
def gene_variants_subset_api(gene_id):
    try:
        intervalset = IntervalSet.from_gene(get_db(), gene_id)
        return _get_variants_subset_response_for_intervalset(intervalset)
    except:_err(); abort(404)
@bp.route('/api/variants/transcript/<transcript_id>', methods=['POST'])
@require_agreement_to_terms_and_store_destination
def transcript_variants_subset_api(transcript_id):
    try:
        intervalset = IntervalSet.from_transcript(get_db(), transcript_id)
        return _get_variants_subset_response_for_intervalset(intervalset)
    except:_err(); abort(404)
@bp.route('/api/variants/region/<chrom>-<start>-<stop>', methods=['POST'])
@require_agreement_to_terms_and_store_destination
def region_variants_subset_api(chrom, start, stop):
    try:
        start,stop = int(start),int(stop); assert stop-start <= MAX_REGION_LENGTH
        intervalset = IntervalSet.from_chrom_start_stop(chrom, start, stop)
        return _get_variants_subset_response_for_intervalset(intervalset)
    except:_err(); abort(404)
def _get_variants_subset_response_for_intervalset(intervalset):
    db = get_db()
    args = json.loads(request.form['args'])
    assert isinstance(args['draw'], int)
    filter_info = json.loads(request.form['filter_info'])
    _log('   '+str(filter_info))
    ret = lookups.get_variants_subset_for_intervalset(
        db, intervalset, args['columns'], args['order'], filter_info, skip=args['start'], length=args['length']
    )
    ret['draw'] = args['draw']
    return jsonify(ret)


@bp.route('/api/coverage/gene/<gene_id>')
@require_agreement_to_terms_and_store_destination
def gene_coverage_api(gene_id):
    try:
        intervalset = IntervalSet.from_gene(get_db(), gene_id)
        return jsonify(get_coverage_handler().get_coverage_for_intervalset(intervalset))
    except:_err(); abort(404)
@bp.route('/api/coverage/transcript/<transcript_id>')
@require_agreement_to_terms_and_store_destination
def transcript_coverage_api(transcript_id):
    try:
        intervalset = IntervalSet.from_transcript(get_db(), transcript_id)
        return jsonify(get_coverage_handler().get_coverage_for_intervalset(intervalset))
    except:_err(); abort(404)
@bp.route('/api/coverage/region/<chrom>-<start>-<stop>')
@require_agreement_to_terms_and_store_destination
def region_coverage_api(chrom, start, stop):
    try:
        start,stop = int(start),int(stop); assert stop-start <= MAX_REGION_LENGTH
        intervalset = IntervalSet.from_chrom_start_stop(chrom, start, stop)
        return jsonify(get_coverage_handler().get_coverage_for_intervalset(intervalset))
    except:_err(); abort(404)


@bp.route('/multi_variant_rsid/<rsid>')
@require_agreement_to_terms_and_store_destination
def multi_variant_rsid_page(rsid):
    db = get_db()
    try:
        _log()
        variants = lookups.get_variants_by_rsid(db, rsid)
        if variants is None or len(variants) == 0:
            return error_page("There are no variants with the rsid '{}'".format(rsid))
        return error_page('There are multiple variants at the location of rsid {}: {}'.format(
            rsid,
            ', '.join('{chrom}-{pos}-{ref}-{alt}'.format(**variant) for variant in variants)))
    except:_err(); abort(404)


@bp.route('/not_found/<query>')
def not_found_page(query):
    return render_template(
        'not_found.html',
        query=query
    )

@bp.route('/error/<message>')
@bp.errorhandler(404)
def error_page(message):
    return render_template(
        'error.html',
        message=message
    ), 404


@bp.route('/download')
@require_agreement_to_terms_and_store_destination
def download_page():
    _log()
    return render_template('download.html')


@bp.route('/download/all')
@require_agreement_to_terms_and_store_destination
def download_full_vcf():
    _log()
    try:
        return make_response(send_file(app.config['DOWNLOAD_ALL_FILEPATH'], as_attachment=True, mimetype='application/gzip'))
    except:_err(); abort(404)


@bp.route('/about')
def about_page():
    _log()
    return render_template('about.html')

@bp.route('/terms')
def terms_page():
    _log()
    return render_template('terms.html')


# OAuth2
google_sign_in = auth.GoogleSignIn(app)

lm = LoginManager(app)
lm.login_view = 'bp.homepage'

class User(UserMixin):
    "A user's id is their email address."
    def __init__(self, username=None, email=None, agreed_to_terms=False, picture=None, enabled_api=False, google_client_id=None, no_newsletters=False):
        self.username = username
        self.email = email
        self.agreed_to_terms = agreed_to_terms
        self.picture = picture
        self.enabled_api = enabled_api
        self.google_client_id = google_client_id
        self.no_newsletters = no_newsletters

    def get_id(self):
        return self.email
    def __str__(self):
        return "<{}>".format(self.email or None)
    def __repr__(self):
        return "<User email={!r} username={!r} terms={!r}>".format(self.email, self.username, self.agreed_to_terms)


def encode_user(user):
    return {'_type': 'User', 'user_id': user.get_id(), 'username': user.username, 'email': user.email, 'agreed_to_terms': user.agreed_to_terms, 'picture': user.picture, 'enabled_api': user.enabled_api, 'google_client_id': user.google_client_id, 'no_newsletters': user.no_newsletters}

def decode_user(document):
    assert document['_type'] == 'User'
    return User(document['username'], document['email'], document['agreed_to_terms'], document.get('picture', None), document.get('enabled_api', False), document.get('google_client_id', None), document.get('no_newsletters', False))

@lm.user_loader
def load_user(id):
    db = get_db()
    document = db.users.find_one({'user_id': id}, projection = {'_id': False})
    if document:
        u = decode_user(document)
    else:
        # This method is supposed to support bad `id`s.
        print('user not found with id [{!r}]'.format(id))
        u = None
    return u

@bp.route('/agree_to_terms')
def agree_to_terms():
    "this route is for when the user has clicked 'I agree to the terms'."
    if not current_user.is_anonymous:
        current_user.agreed_to_terms = True
        db = get_db()
        result = db.users.update_one({"user_id": current_user.get_id()}, {"$set": {"agreed_to_terms": current_user.agreed_to_terms}})
    _log()
    return redirect(url_for('.get_authorized'))

@bp.route('/logout')
def logout():
    _log()
    logout_user()
    return redirect(url_for('.homepage'))

@bp.route('/login_with_google')
def login_with_google():
    "this route is for the login button"
    session['original_destination'] = url_for('.homepage')
    return redirect(url_for('.get_authorized'))

@bp.route('/get_authorized')
def get_authorized():
    "This route tries to be clever and handle lots of situations."
    if current_user.is_anonymous:
        return google_sign_in.authorize()
    elif not current_user.agreed_to_terms:
        return redirect(url_for('.terms_page'))
    else:
        if 'original_destination' in session:
            orig_dest = session['original_destination']
            del session['original_destination'] # We don't want old destinations hanging around.  If this leads to problems with re-opening windows, disable this line.
        else:
            orig_dest = url_for('.homepage')
        return redirect(orig_dest)

@bp.route('/callback/google')
def oauth_callback_google():
    if not current_user.is_anonymous:
        return redirect(url_for('.homepage'))
    try:
        username, email, picture = google_sign_in.callback() # oauth.callback reads request.args.
    except:
        print('Error in google_sign_in.callback():')
        print(traceback.format_exc())
        flash('Something is wrong with authentication.  Please email pjvh@umich.edu')
        return redirect(url_for('.homepage'))
    if email is None:
        # I need a valid email address for my user identification
        flash('Authentication failed by failing to get an email address.  Please email pjvh@umich.edu')
        return redirect(url_for('.homepage'))

    if app.config['EMAIL_WHITELIST']:
        if email.lower() not in app.config['EMAIL_WHITELIST']:
            flash('Your email, {}, is not in the list of allowed emails. If it should be, email pjvh@umich.edu to request permission.'.format(email.lower()))
            return redirect(url_for('.homepage'))

    # Look if the user already exists

    db = get_db()
    document = db.users.find_one({'user_id': email}, projection = {'_id': False})

    if document:
        user = decode_user(document)
        if picture and picture != user.picture:
            result = db.users.update_one({"user_id": user.get_id()}, {"$set": {"picture": picture}})
            user.picture = picture
    else:
        user = User(email=email, username=username or email.split('@')[0], picture=picture)
        db.users.insert(encode_user(user))
    #session['picture'] = None
    # Log in the user, by default remembering them for their next visit
    # unless they log out.
    login_user(user, remember=True)
    return redirect(url_for('.get_authorized'))

@bp.after_request
def apply_caching(response):
    # prevent click-jacking vulnerability identified by BITs
    response.headers["X-Frame-Options"] = "SAMEORIGIN"
    return response

app.register_blueprint(bp, url_prefix = app.config['URL_PREFIX'])

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--host', default='0.0.0.0', help='the hostname to use to access this server')
    parser.add_argument('--port', type=int, default=5000, help='an integer for the accumulator')
    args = parser.parse_args()
    app.run(host=args.host, port=args.port, threaded=True, use_reloader=True)
