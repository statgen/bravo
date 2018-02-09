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
argparser_samples.add_argument('-i', '--in', metavar = 'file', dest = 'inFile', required = True, help = 'Input compressed VCF.')


argparser_bams = subparsers.add_parser('cram', help = 'Generates CRAM file with sequences from heterozygous/homozygous samples.')
argparser_bams.add_argument('-i', '--in', metavar = 'file', dest = 'inFile', required = True, help = 'Input compressed VCF with HET and HOM INFO fields. Multi-allelic variants must be split into bi-allelic entries.')
argparser_bams.add_argument('-c', '--crams', metavar = 'file', dest = 'inCRAMs', required = True, help = 'Input file with sample name and CRAM file path per line.')
argparser_bams.add_argument('-w', '--window', metavar = 'base-pair', dest = 'window', type = int, required = False, default = 100, help = 'Window size around each variant in base-pairs.')
argparser_bams.add_argument('-o', '--out', metavar = 'file', dest = 'outFile', required = True, help = 'Output CRAM. Unsorted.')


Variant = namedtuple('Variant', ['chrom', 'pos', 'ref', 'alt'], verbose = False)
Region = namedtuple('Region', ['variant_idx', 'sample_idx', 'is_hom', 'window'], verbose = False)

variants = []

def process_sample(cram_path, regions, ocram):
    with pysam.AlignmentFile(cram_path, 'rc') as icram:
        for region in regions:
            process_region(icram, region, ocram)


def process_region(icram, region, ocram):
    global variants
    variant = variants[region.variant_idx]
    qnames = dict()
    for read in icram.fetch(variant.chrom, variant.pos - region.window if variant.pos > region.window else 0, variant.pos + region.window):
        qname = qnames.get(read.query_name, None)
        if qname is None:
            qname = str(len(qnames) + 1)
            qnames[read.query_name] = qname
        read.query_name = '{}:{}:{}:{}{}:{}'.format(variant.pos, variant.ref, variant.alt, '0' if region.is_hom else '', region.sample_idx, qname)
        for tag, value in read.get_tags():
            read.set_tag(tag, None)
        ocram.write(read)


if __name__ == '__main__':
    args = argparser.parse_args()

    if args.command == 'samples':
        unique_samples = set()
        with gzip.GzipFile(args.inFile, 'r') as ifile:
            for line in ifile:
                if line.startswith('#'):
                    continue
                fields = line.rstrip().split('\t')
                info = dict(map(lambda x: (x[0], x[1]) if len(x) > 1 else (x[0], None), (x.split('=', 1) for x in fields[7].split(';'))))
                hets = info.get('HET', None)
                homs = info.get('HOM', None)
                if hets is not None:
                    for sample in hets.strip().split(','):
                        unique_samples.add(sample.strip())
                if homs is not None:
                    for sample in homs.strip().split(','):
                        unique_samples.add(sample.strip())
        for sample in unique_samples:
            sys.stdout.write('{}\n'.format(sample))
    elif args.command == 'cram':
        crams = dict()
        header = {'HD': {'SO': 'coordinate', 'VN': '1.3'}, 'SQ': [], 'RG': []}

        with open(args.inCRAMs, 'r') as ifile:
            for line in ifile:
                fields = line.rstrip().split()
                if len(fields) >= 2:
                    sample_id = fields[0].strip()
                    cram_path = fields[1].strip()
                    if os.path.isfile(cram_path):
                        crams[sample_id] = cram_path
                    else:
                        sys.stdout.write('CRAM "{}" for "{}" was not found.\n'.format(cram_path, sample_id))

        if len(crams) == 0:
            sys.exit(0)

        samples = dict()

        with gzip.GzipFile(args.inFile, 'r') as ifile:
            for line in ifile:
                if line.startswith('#'):
                    continue
                fields = line.rstrip().split('\t')
                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                info = dict(map(lambda x: (x[0], x[1]) if len(x) > 1 else (x[0], None), (x.split('=', 1) for x in fields[7].split(';'))))
                hets = info.get('HET', None)
                homs = info.get('HOM', None)
                variant_idx = len(variants)
                variants.append(Variant(chrom, pos, ref, alt))
                sample_idx = 0
                if hets is not None:
                    for sample in (x.strip() for x in hets.strip().split(',')):
                        if sample in crams:
                            if sample not in samples:
                                 samples[sample] = []
                            sample_idx += 1
                            samples[sample].append(Region(variant_idx, sample_idx, False, args.window))
                if homs is not None:
                    for sample in (x.strip() for x in homs.strip().split(',')):
                        if sample in crams:
                            if sample not in samples:
                                samples[sample] = []
                            sample_idx += 1
                            samples[sample].append(Region(variant_idx, sample_idx, True, args.window))

        with pysam.AlignmentFile(crams[next(iter(crams))], 'rc') as icram:
            for sq_line in icram.header['SQ']:
                header['SQ'].append(sq_line)

        i = 0
        with pysam.AlignmentFile(args.outFile, 'wc', header = header) as ocram:
            for sample, sample_variants in samples.iteritems():
                process_sample(crams[sample], sample_variants, ocram)
                i += 1
                if i % 100 == 0:
                    sys.stdout.write('Processed {}/{} sample(s).\n'.format(i, len(samples)))
        sys.stdout.write('Done ({}/{}).\n'.format(i, len(samples)))
