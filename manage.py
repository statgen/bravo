#!/usr/bin/env python2

import sys
import argparse
import imp
import inspect
import warnings
import pymongo
import gzip
import parsing

argparser = argparse.ArgumentParser(description = 'Tool for creating and populating Bravo database.')
argparser_subparsers = argparser.add_subparsers(help = '', dest = 'command')

argparser_gene_models = argparser_subparsers.add_parser('genes', help = 'Creates and populates MongoDB collections for gene models.')
argparser_gene_models.add_argument('-c', '--config', metavar = 'name', required = True, type = str, dest = 'config_class_name', help = 'Bravo configuration class name.')
argparser_gene_models.add_argument('-t', '--canonical-transcripts', metavar = 'file', required = True, type = str, dest = 'canonical_transcripts_file', help = 'File (compressed using Gzip) with a list of canonical transcripts. Must have two columns without a header. First column stores Ensembl gene ID, second column stores Ensembl transcript ID.')
argparser_gene_models.add_argument('-m', '--omim', metavar = 'file', required = True, type = str, dest = 'omim_file', help = 'File (compressed using Gzip) with genes descriptions from OMIM. Required columns separated by tab: Gene stable ID, Transcript stable ID, MIM gene accession, MIM gene description.')
argparser_gene_models.add_argument('-f', '--dbnsfp', metavar = 'file', required = True, type = str, dest = 'genenames_file', help = 'File (compressed using Gzip) with gene names from HGNC. Required columns separated by tab: symbol, name, alias_symbol, prev_name, ensembl_gene_id.')
argparser_gene_models.add_argument('-g', '--gencode', metavar = 'file', required = True, type = str, dest = 'gencode_file', help = 'File from GENCODE in compressed GTF format.')


def load_config(name):
    """Loads Bravo configuration class.
    
    Arguments:
    name -- string in the format '<module_name>.<class_name>'.
    """
    f, path, desc = imp.find_module(name.strip().split('.')[0])
    m = imp.load_module(name, f, path, desc)
    config = dict()
    for member, value in inspect.getmembers(getattr(m, name.strip().split('.')[1])):
        config[member] = value
    return config


def load_gene_models(db, canonical_transcripts_file, omim_file, genenames_file, gencode_file):
    """Creates and populates the following MongoDB collections: genes, transcripts, exons.
    
    Arguments:
    db -- connections to MongoDB database.
    canonical_transcripts_file -- file with a list of canonical transcripts. No header. Two columns: Ensebl gene ID, Ensembl transcript ID.
    omim_file -- file with genes descriptions from OMIM. Required columns separated by tab: Gene stable ID, Transcript stable ID, MIM gene accession, MIM gene description.
    genenames_file -- file with gene names from HGNC. Required columns separated by tab: symbol, name, alias_symbol, prev_name, ensembl_gene_id.
    gencode_file -- file from GENCODE in compressed GTF format. 
    """
    db.genes.drop()
    db.transcripts.drop()
    db.exons.drop()

    canonical_transcripts = dict()
    with gzip.GzipFile(canonical_transcripts_file, 'r') as ifile:
        for gene, transcript in parsing.get_canonical_transcripts(ifile):
            canonical_transcripts[gene] = transcript

    omim_annotations = dict()
    with gzip.GzipFile(omim_file, 'r') as ifile:
        for gene, transcript, accession, description in parsing.get_omim_associations(ifile):
            omim_annotations[gene] = (accession, description) # TODO: what about transcript?

    genenames = dict()
    with gzip.GzipFile(genenames_file, 'r') as ifile:
        for gene in parsing.get_genenames(ifile):
            genenames[gene['ensembl_gene']] = (gene['gene_full_name'], gene['gene_other_names'])

    with gzip.GzipFile(gencode_file, 'r') as ifile:
        for gene in parsing.get_regions_from_gencode_gtf(ifile, {'gene'}):
            gene_id = gene['gene_id']
            if gene_id in canonical_transcripts:
                gene['canonical_transcript'] = canonical_transcripts[gene_id]
            if gene_id in omim_annotations:
                gene['omim_accession'] = omim_annotations[gene_id][0]
                gene['omim_description'] = omim_annotations[gene_id][1]
            if gene_id in genenames:
                gene['full_gene_name'] = genenames[gene_id][0]
                gene['other_names'] = genenames[gene_id][1]
            db.genes.insert_one(gene)
    db.genes.ensure_index('gene_id')
    db.genes.ensure_index('gene_name')
    db.genes.ensure_index('other_names')
    db.genes.ensure_index('xstart')
    db.genes.ensure_index('xstop')
    sys.stdout.write('Inserted {} gene(s).\n'.format(db.genes.count()))
  
    with gzip.GzipFile(gencode_file, 'r') as ifile:
        db.transcripts.insert_many(transcript for transcript in parsing.get_regions_from_gencode_gtf(ifile, {'transcript'}))
    db.transcripts.ensure_index('transcript_id')
    db.transcripts.ensure_index('gene_id')
    sys.stdout.write('Inserted {} transcript(s).\n'.format(db.transcripts.count()))

    with gzip.GzipFile(gencode_file, 'r') as ifile:
        db.exons.insert_many(exon for exon in parsing.get_regions_from_gencode_gtf(ifile, {'exon', 'CDS', 'UTR'}))
    db.exons.ensure_index('exon_id')
    db.exons.ensure_index('transcript_id')
    db.exons.ensure_index('gene_id')
    sys.stdout.write('Inserted {} exon(s).\n'.format(db.exons.count()))


if __name__ == '__main__':
    args = argparser.parse_args()
   
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        config = load_config('flask_config.{}'.format(args.config_class_name))

    mongo_host = config['MONGO']['host']
    mongo_port = config['MONGO']['port']
    mongo_db_name = config['MONGO']['name']

    mongo = pymongo.MongoClient(host = mongo_host, port = mongo_port, connect = False)
    db = mongo[mongo_db_name]

    if args.command == 'genes':
        sys.stdout.write('Start loading genes to {} database.\n'.format(mongo_db_name))
        load_gene_models(db, args.canonical_transcripts_file, args.omim_file, args.genenames_file, args.gencode_file)
        sys.stdout.write('Done loading genes to {} database.\n'.format(mongo_db_name))
    else:
        raise Exception('Command {} is not supported.'.format(args.command))


'''
import sys
from flask_script import Manager
import exac


manager = Manager(exac.app)


@manager.command
def load_variants_file():
    exac.load_variants_file()

@manager.command
def load_dbsnp_file():
    exac.load_dbsnp_file()



@manager.command
def create_users():
    exac.create_users()



@manager.option('-i', '--input', dest = 'metrics_filename', type = str, required = True, help = 'File with the pre-calculated metrics. Every metric must be stored on a separate line in JSON format.')
def load_metrics(metrics_filename):
    "Loads pre-calculated metrics across all variants into 'metrics' Mongo DB collection."
    if not metrics_filename.strip():
        sys.exit("Metrics file name must be a non-empty string.")
    exac.load_metrics(metrics_filename)


@manager.option('-v', '--vcf', dest = 'vcfs', type = str, nargs = '+', required = True, help = 'Input VCF name(s).')
def load_percentiles(vcfs):
    "Loads percentiles for each variant from INFO field in the provided VCF.\
     Percentiles in the INFO field must have '_P' suffix and store two comma separated values: lower bound and upper bound."
    if not all(x.strip() for x in vcfs):
        sys.exit("VCF file name(s) must be a non-empty string(s).")
    exac.load_percentiles(vcfs)


@manager.option('-c', '--collection', dest = 'collection', type = str, required = True, help = 'Destination Mongo collection name.') 
@manager.option('-v', '--vcf', dest = 'vcfs', type = str, nargs = '+', required = True, help = 'Input VCF name(s).')
def load_custom_variants_file(collection, vcfs):
    "Loads variants from the specified VCF file(s) into a new Mongo collection.\
     Useful when there is a need to serve multiple different variants sets (e.g. after subsetting samples) through the API."

    if not collection.strip():
        sys.exit("Collection name must be a non-empty string.")
    if not all(x.strip() for x in vcfs):
        sys.exit("VCF file name(s) must be a non-empty string(s).")
    exac.load_custom_variants_file(collection, vcfs)


@manager.option('-c', '--collection', dest = 'collection', type = str, required = False, help = 'Mongo collection name to store paths to cached BAM/CRAM files.')
def create_sequence_cache(collection):
    "Creates Mongo collection with unique index to store paths to cached BAM\CRAM files for the IGV browser.\
     Important: Mongo will not do any cleaning if cache becomes too large."
    
    if collection is not None and not collection.strip():
        sys.exit("Collection name must be a non-empty string.")
    exac.create_sequence_cache(collection)

@manager.option('-w', '--whitelist', dest = 'whitelist_file', type = str, required = True, help = 'Emails whitelist file. One email per line.')
def load_whitelist(whitelist_file):
    "Creates whitelist Monogo collection and populates it with provided emails."
    if not whitelist_file.strip():
        sys.exit("Whitelist file name must be a non-empty string.")
    exac.load_whitelist(whitelist_file)

if __name__ == "__main__":
    manager.run()
'''
