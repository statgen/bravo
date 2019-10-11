#!/usr/bin/env python3

import sys
import os
import gzip
import argparse
import pysam
from collections import namedtuple
from random import shuffle

argparser = argparse.ArgumentParser(description = '')
subparsers = argparser.add_subparsers(help = '', dest = 'command')

argparser_samples = subparsers.add_parser('samples', help = 'Extracts all unique sample IDs from the HET and HOM INFO fields.')
argparser_samples.add_argument('-i', '--in', metavar = 'file', dest = 'inFile', required = True, help = 'Input compressed VCF. Multi-allelic variants must be split into bi-allelic entries.')


argparser_bams = subparsers.add_parser('cram', help = 'Generates CRAM file with sequences from heterozygous/homozygous samples.')
argparser_bams.add_argument('-i', '--in', metavar = 'file', dest = 'inFile', required = True, help = 'Input compressed VCF with HET and HOM INFO fields. Multi-allelic variants must be split into bi-allelic entries.')
argparser_bams.add_argument('-r', '--reference', metavar = 'file', dest = 'inReference', required = True, help = 'Reference FASTA file.')
argparser_bams.add_argument('-c', '--crams', metavar = 'file', dest = 'inCRAMs', required = True, help = 'Input file with a sample name and the corresponding CRAM file path per line.')
argparser_bams.add_argument('-w', '--window', metavar = 'base-pair', dest = 'window', type = int, required = False, default = 100, help = 'Window size around each variant in base-pairs.')
argparser_bams.add_argument('-o', '--out', metavar = 'file', dest = 'outFile', required = True, help = 'Output CRAM. Unsorted.')


Variant = namedtuple('Variant', ['chrom', 'pos', 'ref', 'alt'])
Region = namedtuple('Region', ['variant_idx', 'sample_idx', 'is_hom'])


variants = []
window_bp = None


def process_sample(cram_path, reference_path, regions, ocram):
    with pysam.AlignmentFile(cram_path, 'rc', reference_filename = reference_path) as icram:
        for region in regions:
            process_region(icram, region, ocram)


def process_region(icram, region, ocram):
    global variants
    variant = variants[region.variant_idx]
    qnames = dict()
    for read in icram.fetch(variant.chrom, max(variant.pos - window_bp, 0), variant.pos + window_bp):
        qname = qnames.get(read.query_name, None)
        if qname is None:
            qname = str(len(qnames) + 1)
            qnames[read.query_name] = qname
        read.query_name = f'{variant.pos}:{variant.ref}:{variant.alt}:{"0" if region.is_hom else ""}{region.sample_idx}:{qname}'
        for tag, value in read.get_tags():
            read.set_tag(tag, None)
        ocram.write(read)


if __name__ == '__main__':
    args = argparser.parse_args()
    window_bp = args.window
    if args.command == 'samples':
        unique_samples = set()
        with pysam.VariantFile(args.inFile, 'r')as ifile:
            for record in ifile:
                for sample in record.info.get('HET', []):
                    unique_samples.add(sample)
                for sample in record.info.get('HOM', []):
                    unique_samples.add(sample)
        for sample in unique_samples:
            sys.stdout.write('{}\n'.format(sample))
    elif args.command == 'cram':
        crams = dict()
        with open(args.inCRAMs, 'r') as ifile:
            for line in ifile:
                fields = line.rstrip().split()
                if len(fields) >= 2:
                    sample_id = fields[0].strip()
                    cram_path = fields[1].strip()
                    crams[sample_id] = cram_path
        if len(crams) == 0:
            sys.exit(0)
        samples = dict()
        max_het = 0
        max_hom = 0
        start = sys.maxsize
        stop = 0
        chrom = set()
        with pysam.VariantFile(args.inFile, 'r') as ifile:
           for record in ifile:
                variant_idx = len(variants)
                variants.append(Variant(record.chrom, record.pos, record.ref, record.alts[0]))
                start = min(start, record.pos)
                stop = max(stop, record.pos)
                chrom.add(record.chrom)
                if len(chrom) > 1:
                    raise Exception(f'Only single chromosome per file is allowed.')
                for idx, sample in enumerate(record.info.get('HET', []), 1):
                    if sample not in crams:
                        raise Exception(f'No CRAM/BAM file for {sample}')
                    samples.setdefault(sample, []).append(Region(variant_idx, idx, False))
                    max_het = max(max_het, idx)
                for idx, sample in enumerate(record.info.get('HOM', []), 1):
                    if sample not in crams:
                        raise Exception(f'No CRAM/BAM file for {sample}')
                    samples.setdefault(sample, []).append(Region(variant_idx, idx, True))
                    max_hom = max(max_hom, idx)

        header = { 'HD': {'SO': 'unsorted', 'VN': '1.6'},
                   'SQ': [],
                   'RG': [],
                   'CO': [ f'MAX_HOM={max_hom};MAX_HET={max_het}', f'REGION={next(iter(chrom))}:{start}-{stop}' ]
                }
        with pysam.AlignmentFile(crams[next(iter(crams))], 'rc') as icram:
            for sq_line in icram.header['SQ']:
                header['SQ'].append(sq_line)
        with pysam.AlignmentFile(args.outFile, 'wc', reference_filename = args.inReference, header = header) as ocram:
            for i, (sample, sample_variants) in enumerate(samples.items(), 1):
                process_sample(crams[sample], args.inReference, sample_variants, ocram)
                if i % 100 == 0:
                    sys.stdout.write('Processed {}/{} sample(s).\n'.format(i, len(samples)))
        sys.stdout.write('Done ({}/{}).\n'.format(i, len(samples)))
