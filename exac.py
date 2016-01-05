#!/usr/bin/env python2

import itertools
import json
import os
import pymongo
import pysam
import gzip
from parsing import *
import lookups
import random
from utils import *
from pycoverage import *
import auth

from flask import Flask, Response, request, session, g, redirect, url_for, abort, render_template, flash, jsonify
from flask.ext.compress import Compress
from flask_errormail import mail_on_500
from flask.ext.login import LoginManager, UserMixin, login_user, logout_user, current_user
from collections import defaultdict
from werkzeug.contrib.cache import NullCache # TODO: for production, use FileSystemCache

from multiprocessing import Process
import glob
import traceback
import time
import sys
import functools

ADMINISTRATORS = (
    'pjvh@umich.edu',
)

app = Flask(__name__)
app.config.from_object('flask_config') # contains `SECRET_KEY`.
mail_on_500(app, ADMINISTRATORS)
Compress(app)
app.config['COMPRESS_DEBUG'] = True
#cache = FileSystemCache('cache_dir', default_timeout=60*60*24*31)
cache = NullCache()

#EXAC_FILES_DIRECTORY = '../exac_data/'
EXAC_FILES_DIRECTORY = '/net/inpsyght/mongo/orig_data/'
REGION_LIMIT = 1E5
EXON_PADDING = 50
# Load default config and override config from an environment variable
app.config.update(dict(
    DB_HOST='topmed.sph.umich.edu',
    DB_PORT=27017,
    DB_NAME='topmed_chr22',
    DEBUG=True,
    LOAD_DB_PARALLEL_PROCESSES = 2,  # contigs assigned to threads, so good to make this a factor of 24 (eg. 2,3,4,6,8)
    #SITES_VCFS=glob.glob(os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'ALL.chr22.*.VEP.vcf.gz')),
    SITES_VCFS=glob.glob(os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'topmed_freeze1b_20151105_4317_snps_indels.chr22.sites.anno.vcf.gz')),
    GENCODE_GTF=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'gencode.gtf.gz'),
    CANONICAL_TRANSCRIPT_FILE=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'canonical_transcripts.txt.gz'),
    OMIM_FILE=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'omim_info.txt.gz'),
    BASE_COVERAGE_FILES=glob.glob(os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'coverage', 'Panel.*.coverage.txt.gz')),
    DBNSFP_FILE=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'dbNSFP2.6_gene.gz'),

    # How to get a snp141.txt.bgz file:
    #   wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp141.txt.gz
    #   zcat snp141.txt.gz | cut -f 1-5 | bgzip -c > snp141.txt.bgz
    #   tabix -0 -s 2 -b 3 -e 4 snp141.txt.bgz

    #   wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/database/organism_data/b142_SNPChrPosOnRef_105.bcp.gz
    #   zcat b142_SNPChrPosOnRef_105.bcp.gz | awk '$3 != ""' | perl -pi -e 's/ +/\t/g' | sort -k2,2 -k3,3n | bgzip -c > dbsnp142.txt.bgz
    #   tabix -s 2 -b 3 -e 3 dbsnp142.txt.bgz
    DBSNP_FILE=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'dbsnp142.txt.bgz')
))

# DT: load tabix files for base coverage
FULL_COVERAGE_FILES = {
   '1' : '/net/topmed/working/mongo/imported_data/coverage/full_json/1.topmed_freeze1c_20151119_4317.coverage.json.gz', 
   '2' : '/net/topmed/working/mongo/imported_data/coverage/full_json/2.topmed_freeze1c_20151119_4317.coverage.json.gz', 
   '3' : '/net/topmed/working/mongo/imported_data/coverage/full_json/3.topmed_freeze1c_20151119_4317.coverage.json.gz',
   '4' : '/net/topmed/working/mongo/imported_data/coverage/full_json/4.topmed_freeze1c_20151119_4317.coverage.json.gz', 
   '5' : '/net/topmed/working/mongo/imported_data/coverage/full_json/5.topmed_freeze1c_20151119_4317.coverage.json.gz', 
   '6' : '/net/topmed/working/mongo/imported_data/coverage/full_json/6.topmed_freeze1c_20151119_4317.coverage.json.gz', 
   '7' : '/net/topmed/working/mongo/imported_data/coverage/full_json/7.topmed_freeze1c_20151119_4317.coverage.json.gz', 
   '8' : '/net/topmed/working/mongo/imported_data/coverage/full_json/8.topmed_freeze1c_20151119_4317.coverage.json.gz', 
   '9' : '/net/topmed/working/mongo/imported_data/coverage/full_json/9.topmed_freeze1c_20151119_4317.coverage.json.gz', 
   '10' : '/net/topmed/working/mongo/imported_data/coverage/full_json/10.topmed_freeze1c_20151119_4317.coverage.json.gz',
   '11' : '/net/topmed/working/mongo/imported_data/coverage/full_json/11.topmed_freeze1c_20151119_4317.coverage.json.gz',
   '12' : '/net/topmed/working/mongo/imported_data/coverage/full_json/12.topmed_freeze1c_20151119_4317.coverage.json.gz',
   '13' : '/net/topmed/working/mongo/imported_data/coverage/full_json/13.topmed_freeze1c_20151119_4317.coverage.json.gz',
   '14' : '/net/topmed/working/mongo/imported_data/coverage/full_json/14.topmed_freeze1c_20151119_4317.coverage.json.gz',
   '15' : '/net/topmed/working/mongo/imported_data/coverage/full_json/15.topmed_freeze1c_20151119_4317.coverage.json.gz',
   '16' : '/net/topmed/working/mongo/imported_data/coverage/full_json/16.topmed_freeze1c_20151119_4317.coverage.json.gz',
   '17' : '/net/topmed/working/mongo/imported_data/coverage/full_json/17.topmed_freeze1c_20151119_4317.coverage.json.gz',
   '18' : '/net/topmed/working/mongo/imported_data/coverage/full_json/18.topmed_freeze1c_20151119_4317.coverage.json.gz',
   '19' : '/net/topmed/working/mongo/imported_data/coverage/full_json/19.topmed_freeze1c_20151119_4317.coverage.json.gz',
   '20' : '/net/topmed/working/mongo/imported_data/coverage/full_json/20.topmed_freeze1c_20151119_4317.coverage.json.gz',
   '21' : '/net/topmed/working/mongo/imported_data/coverage/full_json/21.topmed_freeze1c_20151119_4317.coverage.json.gz',
   '22' : '/net/topmed/working/mongo/imported_data/coverage/full_json/22.topmed_freeze1c_20151119_4317.coverage.json.gz'
}

BINNED_COVERAGE_FILES = {
   '1' : '/net/topmed/working/mongo/imported_data/coverage/bins_json/1.topmed_freeze1c_20151119_4317.coverage.bins.json.gz', 
   '2' : '/net/topmed/working/mongo/imported_data/coverage/bins_json/2.topmed_freeze1c_20151119_4317.coverage.bins.json.gz', 
   '3' : '/net/topmed/working/mongo/imported_data/coverage/bins_json/3.topmed_freeze1c_20151119_4317.coverage.bins.json.gz',
   '4' : '/net/topmed/working/mongo/imported_data/coverage/bins_json/4.topmed_freeze1c_20151119_4317.coverage.bins.json.gz', 
   '5' : '/net/topmed/working/mongo/imported_data/coverage/bins_json/5.topmed_freeze1c_20151119_4317.coverage.bins.json.gz', 
   '6' : '/net/topmed/working/mongo/imported_data/coverage/bins_json/6.topmed_freeze1c_20151119_4317.coverage.bins.json.gz', 
   '7' : '/net/topmed/working/mongo/imported_data/coverage/bins_json/7.topmed_freeze1c_20151119_4317.coverage.bins.json.gz', 
   '8' : '/net/topmed/working/mongo/imported_data/coverage/bins_json/8.topmed_freeze1c_20151119_4317.coverage.bins.json.gz', 
   '9' : '/net/topmed/working/mongo/imported_data/coverage/bins_json/9.topmed_freeze1c_20151119_4317.coverage.bins.json.gz', 
   '10' : '/net/topmed/working/mongo/imported_data/coverage/bins_json/10.topmed_freeze1c_20151119_4317.coverage.bins.json.gz',
   '11' : '/net/topmed/working/mongo/imported_data/coverage/bins_json/11.topmed_freeze1c_20151119_4317.coverage.bins.json.gz',
   '12' : '/net/topmed/working/mongo/imported_data/coverage/bins_json/12.topmed_freeze1c_20151119_4317.coverage.bins.json.gz',
   '13' : '/net/topmed/working/mongo/imported_data/coverage/bins_json/13.topmed_freeze1c_20151119_4317.coverage.bins.json.gz',
   '14' : '/net/topmed/working/mongo/imported_data/coverage/bins_json/14.topmed_freeze1c_20151119_4317.coverage.bins.json.gz',
   '15' : '/net/topmed/working/mongo/imported_data/coverage/bins_json/15.topmed_freeze1c_20151119_4317.coverage.bins.json.gz',
   '16' : '/net/topmed/working/mongo/imported_data/coverage/bins_json/16.topmed_freeze1c_20151119_4317.coverage.bins.json.gz',
   '17' : '/net/topmed/working/mongo/imported_data/coverage/bins_json/17.topmed_freeze1c_20151119_4317.coverage.bins.json.gz',
   '18' : '/net/topmed/working/mongo/imported_data/coverage/bins_json/18.topmed_freeze1c_20151119_4317.coverage.bins.json.gz',
   '19' : '/net/topmed/working/mongo/imported_data/coverage/bins_json/19.topmed_freeze1c_20151119_4317.coverage.bins.json.gz',
   '20' : '/net/topmed/working/mongo/imported_data/coverage/bins_json/20.topmed_freeze1c_20151119_4317.coverage.bins.json.gz',
   '21' : '/net/topmed/working/mongo/imported_data/coverage/bins_json/21.topmed_freeze1c_20151119_4317.coverage.bins.json.gz',
   '22' : '/net/topmed/working/mongo/imported_data/coverage/bins_json/22.topmed_freeze1c_20151119_4317.coverage.bins.json.gz'
}

coverages = CoverageCollection()

for contig, file in FULL_COVERAGE_FILES.iteritems():
    coverages.setTabixPath(0, 50000, contig, file)
#    coverages.setTabixPath(0, sys.maxint, contig, file)

for contig, file in BINNED_COVERAGE_FILES.iteritems():
    coverages.setTabixPath(50000, sys.maxint, contig, file)

coverages.openAll()

def connect_db():
    """
    Connects to the specific database.
    """
    client = pymongo.MongoClient(host=app.config['DB_HOST'], port=app.config['DB_PORT'])
    return client[app.config['DB_NAME']]


def parse_tabix_file_subset(tabix_filenames, subset_i, subset_n, record_parser):
    """
    Returns a generator of parsed record objects (as returned by record_parser) for the i'th out n subset of records
    across all the given tabix_file(s). The records are split by files and contigs within files, with 1/n of all contigs
    from all files being assigned to this the i'th subset.

    Args:
        tabix_filenames: a list of one or more tabix-indexed files. These will be opened using pysam.Tabixfile
        subset_i: zero-based number
        subset_n: total number of subsets
        record_parser: a function that takes a file-like object and returns a generator of parsed records
    """
    start_time = time.time()
    open_tabix_files = [pysam.Tabixfile(tabix_filename) for tabix_filename in tabix_filenames]
    tabix_file_contig_pairs = [(tabix_file, contig) for tabix_file in open_tabix_files for contig in tabix_file.contigs]
    tabix_file_contig_subset = tabix_file_contig_pairs[subset_i : : subset_n]  # get every n'th tabix_file/contig pair
    short_filenames = ", ".join(map(os.path.basename, tabix_filenames))
    num_file_contig_pairs = len(tabix_file_contig_subset)
    print(("Loading subset %(subset_i)s of %(subset_n)s total: %(num_file_contig_pairs)s contigs from "
           "%(short_filenames)s") % locals())
    counter = 0
    for tabix_file, contig in tabix_file_contig_subset:
        header_iterator = tabix_file.header
        records_iterator = tabix_file.fetch(contig, 0, 10**9, multiple_iterators=True)
        for parsed_record in record_parser(itertools.chain(header_iterator, records_iterator)):
            counter += 1
            yield parsed_record

            if counter % 100000 == 0:
                seconds_elapsed = int(time.time()-start_time)
                print(("Loaded %(counter)s records from subset %(subset_i)s of %(subset_n)s from %(short_filenames)s "
                       "(%(seconds_elapsed)s seconds)") % locals())

    print("Finished loading subset %(subset_i)s from  %(short_filenames)s (%(counter)s records)" % locals())


def load_base_coverage():
    def load_coverage(coverage_files, i, n, db):
        coverage_generator = parse_tabix_file_subset(coverage_files, i, n, get_base_coverage_from_file)
        try:
            db.base_coverage.insert(coverage_generator, w=0)
        except pymongo.errors.InvalidOperation:
            pass  # handle error when coverage_generator is empty

    db = get_db()
    db.base_coverage.drop()
    print("Dropped db.base_coverage")
    # load coverage first; variant info will depend on coverage
    db.base_coverage.ensure_index('xpos')

    procs = []
    coverage_files = app.config['BASE_COVERAGE_FILES']
    num_procs = app.config['LOAD_DB_PARALLEL_PROCESSES']
    random.shuffle(app.config['BASE_COVERAGE_FILES'])
    for i in range(num_procs):
        p = Process(target=load_coverage, args=(coverage_files, i, num_procs, db))
        p.start()
        procs.append(p)
    return procs

    #print 'Done loading coverage. Took %s seconds' % int(time.time() - start_time)


def load_variants_file():
    def load_variants(sites_file, i, n, db):
        variants_generator = parse_tabix_file_subset([sites_file], i, n, get_variants_from_sites_vcf)
        try:
            db.variants.insert(variants_generator, w=0)
        except pymongo.errors.InvalidOperation:
            pass  # handle error when variant_generator is empty

    db = get_db()
    db.variants.drop()
    print("Dropped db.variants")

    # grab variants from sites VCF
    db.variants.ensure_index('xpos')
    db.variants.ensure_index('xstart')
    db.variants.ensure_index('xstop')
    db.variants.ensure_index('rsids')
    db.variants.ensure_index('genes')
    db.variants.ensure_index('transcripts')

    sites_vcfs = app.config['SITES_VCFS']
    if len(sites_vcfs) > 1:
        raise Exception("More than one sites vcf file found: %s" % sites_vcfs)

    procs = []
    num_procs = app.config['LOAD_DB_PARALLEL_PROCESSES']
    for i in range(num_procs):
        p = Process(target=load_variants, args=(sites_vcfs[0], i, num_procs, db))
        p.start()
        procs.append(p)
    return procs

    #print 'Done loading variants. Took %s seconds' % int(time.time() - start_time)


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

    return []


def load_dbsnp_file():
    db = get_db()

    def load_dbsnp(dbsnp_file, i, n, db):
        if os.path.isfile(dbsnp_file + ".tbi"):
            dbsnp_record_generator = parse_tabix_file_subset([dbsnp_file], i, n, get_snp_from_dbsnp_file)
            try:
                db.dbsnp.insert(dbsnp_record_generator, w=0)
            except pymongo.errors.InvalidOperation:
                pass  # handle error when coverage_generator is empty

        else:
            with gzip.open(dbsnp_file) as f:
                db.dbsnp.insert((snp for snp in get_snp_from_dbsnp_file(f)), w=0)

    db.dbsnp.drop()
    db.dbsnp.ensure_index('rsid')
    db.dbsnp.ensure_index('xpos')
    start_time = time.time()
    dbsnp_file = app.config['DBSNP_FILE']

    print "Loading dbsnp from %s" % dbsnp_file
    if os.path.isfile(dbsnp_file + ".tbi"):
        num_procs = app.config['LOAD_DB_PARALLEL_PROCESSES']
    else:
        # see if non-tabixed .gz version exists
        if os.path.isfile(dbsnp_file):
            print(("WARNING: %(dbsnp_file)s.tbi index file not found. Will use single thread to load dbsnp."
                "To create a tabix-indexed dbsnp file based on UCSC dbsnp, do: \n"
                "   wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp141.txt.gz \n"
                "   gzcat snp141.txt.gz | cut -f 1-5 | bgzip -c > snp141.txt.bgz \n"
                "   tabix -0 -s 2 -b 3 -e 4 snp141.txt.bgz") % locals())
            num_procs = 1
        else:
            raise Exception("dbsnp file %s(dbsnp_file)s not found." % locals())

    procs = []
    for i in range(num_procs):
        p = Process(target=load_dbsnp, args=(dbsnp_file, i, num_procs, db))
        p.start()
        procs.append(p)

    return procs
    #print 'Done loading dbSNP. Took %s seconds' % int(time.time() - start_time)

    #start_time = time.time()
    #db.dbsnp.ensure_index('rsid')
    #print 'Done indexing dbSNP table. Took %s seconds' % int(time.time() - start_time)


def load_db():
    """
    Load the database
    """
    # Initialize database
    # Don't need to explicitly create tables with mongo, just indices
    confirm = raw_input('This will drop the database and reload. Are you sure you want to continue? [no] ')
    if not confirm.startswith('y'):
        print('Exiting...')
        sys.exit(1)
    all_procs = []
    for load_function in [load_variants_file, load_dbsnp_file, load_base_coverage, load_gene_models]:
        procs = load_function()
        all_procs.extend(procs)
        print("Started %s processes to run %s" % (len(procs), load_function.__name__))

    [p.join() for p in all_procs]
    print('Done! Creating cache...')
    create_cache()
    print('Done!')


def create_cache():
    """
    This is essentially a compile step that generates all cached resources.
    Creates files like autocomplete_entries.txt
    Should be run on every redeploy.
    """
    # create autocomplete_entries.txt
    autocomplete_strings = []
    print >> sys.stderr, "Getting gene names..."
    for gene in get_db().genes.find():
        autocomplete_strings.append(gene['gene_name'])
        if 'other_names' in gene:
            autocomplete_strings.extend(gene['other_names'])
    print >> sys.stderr, "Done! Writing..."
    f = open(os.path.join(os.path.dirname(__file__), 'autocomplete_strings.txt'), 'w')
    for s in sorted(autocomplete_strings):
        f.write(s+'\n')
    f.close()
    print >> sys.stderr, "Done!"


def precalculate_variant_consqequence_category():
    db = get_db()
    print 'Reading %s variants' % db.variants.count()
    for variant in db.variants.find(projection=['_id', 'vep_annotations']):
        add_consequence_to_variant(variant)
        
    db.variants.ensure_index('consequence_category')

def precalculate_metrics():
    import numpy
    db = get_db()
    print 'Reading %s variants...' % db.variants.count()
    metrics = defaultdict(list)
    binned_metrics = defaultdict(list)
    progress = 0
    start_time = time.time()
    for variant in db.variants.find(projection=['quality_metrics', 'site_quality', 'allele_num', 'allele_count']):
        for metric, value in variant['quality_metrics'].iteritems():
            metrics[metric].append(float(value))
        qual = float(variant['site_quality'])
        metrics['site_quality'].append(qual)
        if variant['allele_num'] == 0: continue
        if variant['allele_count'] == 1:
            binned_metrics['singleton'].append(qual)
        elif variant['allele_count'] == 2:
            binned_metrics['doubleton'].append(qual)
        else:
            for af in AF_BUCKETS:
                if float(variant['allele_count'])/variant['allele_num'] < af:
                    binned_metrics[af].append(qual)
                    break
        progress += 1
        if not progress % 100000:
            print 'Read %s variants. Took %s seconds' % (progress, int(time.time() - start_time))
    print 'Done reading variants. Dropping metrics database... '
    db.metrics.drop()
    print 'Dropped metrics database. Calculating metrics...'
    for metric in metrics:
        bin_range = None
        data = map(numpy.log, metrics[metric]) if metric == 'DP' else metrics[metric]
        if metric == 'FS':
            bin_range = (0, 20)
        elif metric == 'VQSLOD':
            bin_range = (-20, 20)
        elif metric == 'InbreedingCoeff':
            bin_range = (0, 1)
        if bin_range is not None:
            data = [x if (x > bin_range[0]) else bin_range[0] for x in data]
            data = [x if (x < bin_range[1]) else bin_range[1] for x in data]
        hist = numpy.histogram(data, bins=40, range=bin_range)
        edges = hist[1]
        # mids = [(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)]
        lefts = [edges[i] for i in range(len(edges)-1)]
        db.metrics.insert({
            'metric': metric,
            'mids': lefts,
            'hist': list(hist[0])
        })
    for metric in binned_metrics:
        hist = numpy.histogram(map(numpy.log, binned_metrics[metric]), bins=40)
        edges = hist[1]
        mids = [(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)]
        db.metrics.insert({
            'metric': 'binned_%s' % metric,
            'mids': mids,
            'hist': list(hist[0])
        })
    db.metrics.ensure_index('metric')
    print 'Done pre-calculating metrics!'


def get_db():
    """
    Opens a new database connection if there is none yet for the
    current application context.
    """
    if not hasattr(g, 'db_conn'):
        g.db_conn = connect_db()
    return g.db_conn


# @app.teardown_appcontext
# def close_db(error):
#     """Closes the database again at the end of the request."""
#     if hasattr(g, 'db_conn'):
#         g.db_conn.close()

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
            return redirect(url_for('get_authorized'))
        return func(*args, **kwargs)
    return decorated_view


@app.route('/')
def homepage():
    cache_key = 't-homepage'
    t = cache.get(cache_key)
    if t is None:
        t = render_template('homepage.html')
        cache.set(cache_key, t)
    return t


@app.route('/autocomplete/<query>')
def awesome_autocomplete(query):
    if not hasattr(g, 'autocomplete_strings'):
        g.autocomplete_strings = [s.strip() for s in open(os.path.join(os.path.dirname(__file__), 'autocomplete_strings.txt'))]
    suggestions = lookups.get_awesomebar_suggestions(g, query)
    return Response(json.dumps([{'value': s} for s in suggestions]),  mimetype='application/json')


@app.route('/awesome')
def awesome():
    db = get_db()
    query = request.args.get('query')
    datatype, identifier = lookups.get_awesomebar_result(db, query)

    print "Searched for %s: %s" % (datatype, identifier)
    if datatype == 'gene':
        return redirect('/gene/{}'.format(identifier))
    elif datatype == 'transcript':
        return redirect('/transcript/{}'.format(identifier))
    elif datatype == 'variant':
        return redirect('/variant/{}'.format(identifier))
    elif datatype == 'region':
        return redirect('/region/{}'.format(identifier))
    elif datatype == 'dbsnp_variant_set':
        return redirect('/dbsnp/{}'.format(identifier))
    elif datatype == 'not_found':
        return redirect('/not_found/{}'.format(identifier))
    else:
        raise Exception


@app.route('/variant/<variant_str>')
@require_agreement_to_terms_and_store_destination
def variant_page(variant_str):
    db = get_db()
    try:
        chrom, pos, ref, alt = variant_str.split('-')
        pos = int(pos)
        # pos, ref, alt = get_minimal_representation(pos, ref, alt)
        xpos = get_xpos(chrom, pos)
        variant = lookups.get_variant(db, xpos, ref, alt)

        if variant is None:
            variant = {
                'chrom': chrom,
                'pos': pos,
                'xpos': xpos,
                'ref': ref,
                'alt': alt
            }
        consequences = None
        ordered_csqs = None
        if 'vep_annotations' in variant:
            variant['vep_annotations'] = order_vep_by_csq(variant['vep_annotations'])  # Adds major_consequence
            ordered_csqs = [x['major_consequence'] for x in variant['vep_annotations']]
            ordered_csqs = reduce(lambda x, y: ','.join([x, y]) if y not in x else x, ordered_csqs, '').split(',') # Close but not quite there
            consequences = defaultdict(lambda: defaultdict(list))
            for annotation in variant['vep_annotations']:
                annotation['HGVS'] = get_proper_hgvs(annotation)
                consequences[annotation['major_consequence']][annotation['Gene']].append(annotation)
        # DT: get coverage from tabix
        base_coverage = lookups.get_coverage_for_bases(coverages, xpos, xpos + len(ref) - 1)
        #base_coverage = lookups.get_coverage_for_bases(db, xpos, xpos + len(ref) - 1)
        any_covered = any([x['has_coverage'] for x in base_coverage])
        metrics = lookups.get_metrics(db, variant)

        print 'Rendering variant: %s' % variant_str
        return render_template(
            'variant.html',
            variant=variant,
            base_coverage=base_coverage,
            consequences=consequences,
            any_covered=any_covered,
            ordered_csqs=ordered_csqs,
            metrics=metrics
        )
    except Exception, e:
        print 'Failed on variant:', variant_str, '; Error=', traceback.format_exc()
        abort(404)


@app.route('/gene/<gene_id>')
@require_agreement_to_terms_and_store_destination
def gene_page(gene_id):
    return get_gene_page_content(gene_id)

def get_gene_page_content(gene_id):
    db = get_db()
    try:
        gene = lookups.get_gene(db, gene_id)
        if gene is None:
            abort(404)
        cache_key = 't-gene-{}'.format(gene_id)
        t = cache.get(cache_key)
        print 'Rendering %sgene: %s' % ('' if t is None else 'cached ', gene_id)
        if t is None:
            variants_in_gene = lookups.get_most_important_variants_in_gene(db, gene_id, limit=200)
            num_variants_in_gene = lookups.get_num_variants_in_gene(db, gene_id)
            transcripts_in_gene = lookups.get_transcripts_in_gene(db, gene_id)

            # Get some canonical transcript and corresponding info
            transcript_id = gene['canonical_transcript']
            transcript = lookups.get_transcript(db, transcript_id)
            variants_in_transcript = lookups.get_most_important_variants_in_transcript(db, transcript_id, limit=200)
            # DT: get coverage from tabix
            coverage_stats = lookups.get_coverage_for_transcript(coverages, transcript['xstart'] - EXON_PADDING, transcript['xstop'] + EXON_PADDING) #, num_bins=1000)
            #coverage_stats = lookups.get_coverage_for_transcript(db, transcript['xstart'] - EXON_PADDING, transcript['xstop'] + EXON_PADDING) #, num_bins=1000)
            add_transcript_coordinate_to_variants(db, variants_in_transcript, transcript_id)

            t = render_template(
                'gene.html',
                gene=gene,
                transcript=transcript,
                variants_in_gene=variants_in_gene,
                num_variants_in_gene=num_variants_in_gene,
                variants_in_transcript=variants_in_transcript,
                transcripts_in_gene=transcripts_in_gene,
                coverage_stats=coverage_stats,
                csq_order=csq_order,
            )
            cache.set(cache_key, t)
        return t
    except Exception, e:
        print 'Failed on gene:', gene_id, ';Error=', traceback.format_exc()
        abort(404)


@app.route('/transcript/<transcript_id>')
@require_agreement_to_terms_and_store_destination
def transcript_page(transcript_id):
    db = get_db()
    try:
        transcript = lookups.get_transcript(db, transcript_id)

        cache_key = 't-transcript-{}'.format(transcript_id)
        t = cache.get(cache_key)
        print 'Rendering %stranscript: %s' % ('' if t is None else 'cached ', transcript_id)
        if t is None:

            gene = lookups.get_gene(db, transcript['gene_id'])
            gene['transcripts'] = lookups.get_transcripts_in_gene(db, transcript['gene_id'])
            variants_in_transcript = lookups.get_variants_in_transcript(db, transcript_id)
            # DT: get coverage from tabix
            coverage_stats = lookups.get_coverage_for_transcript(coverages, transcript['xstart'] - EXON_PADDING, transcript['xstop'] + EXON_PADDING)
            #coverage_stats = lookups.get_coverage_for_transcript(db, transcript['xstart'] - EXON_PADDING, transcript['xstop'] + EXON_PADDING)

            add_transcript_coordinate_to_variants(db, variants_in_transcript, transcript_id)

            t = render_template(
                'transcript.html',
                transcript=transcript,
                transcript_json=json.dumps(transcript),
                variants_in_transcript=variants_in_transcript,
                variants_in_transcript_json=json.dumps(variants_in_transcript),
                coverage_stats=coverage_stats,
                coverage_stats_json=json.dumps(coverage_stats),
                gene=gene,
                gene_json=json.dumps(gene),
                csq_order=csq_order,
            )
            
            cache.set(cache_key, t)
        return t
    except Exception, e:
        print 'Failed on transcript:', transcript_id, ';Error=', traceback.format_exc()
        abort(404)

@app.route('/api/variants_in_gene/<gene_id>')
@require_agreement_to_terms_and_store_destination
def variants_gene_api(gene_id):
    # TODO use `cache`
    db = get_db()
    try:
        variants_in_gene = lookups.get_variants_in_gene(db, gene_id)
        return Response(json.dumps(variants_in_gene), mimetype='application/json')
    except Exception as e:
        print 'Failed on gene:', gene_id, ';Error=', traceback.format_exc()
        abort(404)

@app.route('/api/variants_in_transcript/<transcript_id>')
@require_agreement_to_terms_and_store_destination
def variants_transcript_api(transcript_id):
    # TODO use `cache`
    db = get_db()
    try:
        variants_in_transcript = lookups.get_variants_in_transcript(db, transcript_id)
        add_transcript_coordinate_to_variants(db, variants_in_transcript, transcript_id)
        return Response(json.dumps(variants_in_transcript), mimetype='application/json')
    except Exception as e:
        print 'Failed on transcript:', transcript_id, ';Error=', traceback.format_exc()
        abort(404)

@app.route('/api/variants_in_region/<region_id>')
@require_agreement_to_terms_and_store_destination
def variants_region_api(region_id):
    # TODO use `cache`
    db = get_db()
    try:
        chrom, start, stop = region_id.split('-')
        start, stop = int(start), int(stop)
        variants_in_region = lookups.get_variants_in_region(db, chrom, start, stop)
        return Response(json.dumps(variants_in_region), mimetype='application/json')
    except Exception as e:
        print 'Failed on region:', region_id, ';Error=', traceback.format_exc()
        abort(404)

@app.route('/region/<region_id>')
@require_agreement_to_terms_and_store_destination
def region_page(region_id):
    db = get_db()
    try:
        region = region_id.split('-')
        cache_key = 't-region-{}'.format(region_id)
        t = cache.get(cache_key)
        print 'Rendering %sregion: %s' % ('' if t is None else 'cached ', region_id)
        if t is None:
            chrom = region[0]
            start = None
            stop = None
            if len(region) == 3:
                chrom, start, stop = region
                start = int(start)
                stop = int(stop)
            if start is None or stop - start > REGION_LIMIT or stop < start:
                return render_template(
                    'region.html',
                    genes_in_region=None,
                    variants_in_region=None,
                    chrom=chrom,
                    start=start,
                    stop=stop,
                    coverage=None,
                    csq_order=csq_order,
                )
            if start == stop:
                start -= 20
                stop += 20
            genes_in_region = lookups.get_genes_in_region(db, chrom, start, stop)
            variants_in_region = lookups.get_variants_in_region(db, chrom, start, stop)
            xstart = get_xpos(chrom, start)
            xstop = get_xpos(chrom, stop)
            # DT: get coverage from tabix
            coverage_array = lookups.get_coverage_for_bases(coverages, xstart, xstop)
            #coverage_array = lookups.get_coverage_for_bases(db, xstart, xstop)
            t = render_template(
                'region.html',
                genes_in_region=genes_in_region,
                variants_in_region=variants_in_region,
                chrom=chrom,
                start=start,
                stop=stop,
                coverage=coverage_array,
                csq_order=csq_order,
            )
            cache.set(cache_key, t)
        return t
    except Exception, e:
        print 'Failed on region:', region_id, ';Error=', traceback.format_exc()
        abort(404)


@app.route('/dbsnp/<rsid>')
@require_agreement_to_terms_and_store_destination
def dbsnp_page(rsid):
    db = get_db()
    try:
        variants = lookups.get_variants_by_rsid(db, rsid)
        chrom = None
        start = None
        stop = None
        print 'Rendering rsid: %s' % rsid
        return render_template(
            'region.html',
            rsid=rsid,
            variants_in_region=variants,
            chrom=chrom,
            start=start,
            stop=stop,
            coverage=None,
            genes_in_region=None,
            csq_order=csq_order,
        )
    except Exception, e:
        print 'Failed on rsid:', rsid, ';Error=', traceback.format_exc()
        abort(404)


@app.route('/not_found/<query>')
def not_found_page(query):
    return render_template(
        'not_found.html',
        query=query
    )


@app.route('/error/<query>')
@app.errorhandler(404)
def error_page(query):
    return render_template(
        'error.html',
        query=query
    )


@app.route('/about')
def about_page():
    return render_template('about.html')


@app.route('/participants')
def participants_page():
    return render_template('about.html')


@app.route('/terms')
def terms_page():
    return render_template('terms.html')


@app.route('/contact')
def contact_page():
    return render_template('contact.html')


@app.route('/faq')
def faq_page():
    return render_template('faq.html')


@app.route('/text')
def text_page():
    db = get_db()
    query = request.args.get('text')
    datatype, identifier = lookups.get_awesomebar_result(db, query)
    if datatype in ['gene', 'transcript']:
        gene = lookups.get_gene(db, identifier)
        link = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr%(chrom)s%%3A%(start)s-%(stop)s" % gene
        output = '''Searched for %s. Found %s.
%s; Canonical: %s.
%s''' % (query, identifier, gene['full_gene_name'], gene['canonical_transcript'], link)
        output += '' if 'omim_accession' not in gene else '''
In OMIM: %(omim_description)s
http://omim.org/entry/%(omim_accession)s''' % gene
        return output
    elif datatype == 'error' or datatype == 'not_found':
        return "Gene/transcript %s not found" % query
    else:
        return "Search types other than gene transcript not yet supported"

# OAuth2
users = {} # email -> User
google_sign_in = auth.GoogleSignIn(app)

lm = LoginManager(app)
lm.login_view = 'homepage'

class User(UserMixin):
    "A user's id is their email address."
    def __init__(self, username=None, email=None):
        self.username = username
        self.email = email
        self.agreed_to_terms = False
    def get_id(self):
        return self.email
    def __repr__(self):
        return "<User email={!r} username={!r} agreed_to_terms={!r}>".format(self.email, self.username, self.agreed_to_terms)

@lm.user_loader
def load_user(id):
    try:
        u = users[id]
        print('user [{!r}] found with id [{!r}]'.format(u, id))
    except KeyError:
        print('user not found with id [{!r}]'.format(id))
        u = None # This method is supposed to support bad `id`s.
    return u

@app.route('/agreed_to_terms')
def agree_to_terms():
    "this route is for when the user has clicked 'I agree to the terms'."
    if not current_user.is_anonymous:
        current_user.agreed_to_terms = True
    print('User [{!r}] agreed to the terms!'.format(current_user))
    return redirect(url_for('get_authorized'))

@app.route('/logout')
def logout():
    print('logging out user {!r}'.format(current_user))
    logout_user()
    return redirect(url_for('homepage'))

@app.route('/login_with_google')
def login_with_google():
    "this route is for the login button"
    session['original_destination'] = url_for('homepage')
    return redirect(url_for('get_authorized'))

@app.route('/get_authorized')
def get_authorized():
    "This route tries to be clever and handle lots of situations."
    if current_user.is_anonymous:
        return google_sign_in.authorize()
    elif not current_user.agreed_to_terms:
        return redirect(url_for('terms_page'))
    else:
        if 'original_destination' in session:
            orig_dest = session['original_destination']
            del session['original_destination'] # We don't want old destinations hanging around.  If this leads to problems with re-opening windows, disable this line.
        else:
            orig_dest = '/'
        return redirect(orig_dest)

@app.route('/callback/google')
def oauth_callback_google():
    if not current_user.is_anonymous:
        return redirect(url_for('homepage'))
    username, email = google_sign_in.callback() # oauth.callback reads request.args.
    if email is None:
        # I need a valid email address for my user identification
        flash('Authentication failed.')
        return redirect(url_for('homepage'))
    # Look if the user already exists
    try:
        user = users[email]
    except KeyError:
        user = User(email=email, username=username or email.split('@')[0])
        users[email] = user

    # Log in the user, by default remembering them for their next visit
    # unless they log out.
    login_user(user, remember=True)

    return redirect(url_for('get_authorized'))

if __name__ == "__main__":
    app.run(host='localhost', port=7777)
