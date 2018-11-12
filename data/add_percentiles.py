import sys
import gzip
import pysam
import re
from contextlib import closing
import argparse

argparser = argparse.ArgumentParser(description = 'Adds CADD scores to VCF.')
argparser.add_argument('-i', '--in-vcf', metavar = 'file', dest = 'in_VCF', required = True, help = 'Input VCF/BCF file. File may be compressed with gzip/bgzip.')
argparser.add_argument('-p', '--in-percentiles', metavar = 'file', dest = 'in_pctl_files', nargs = '+', required = True, help = 'Input VCF(s) with percentile in INFO field. Must be compressed with bgzip and indexed with tabix. Percentile INFO field must have _PCTL suffix.')
argparser.add_argument('-o', '--out-vcf', metavar = 'file', dest = 'out_VCF', required = True, help = 'Output VCF file. File will be compressed with bgzip.')


class Percentiles(object):
    def __init__(self, VCF):
        self.step_bp = 1000000
        self.VCF = VCF
        self.has_chr_prefix = False
        with pysam.Tabixfile(self.VCF, 'r') as tabix:
            if any(c.startswith('chr') for c in tabix.contigs):
                self.has_chr_prefix = True
        self.fields = []
        with pysam.VariantFile(self.VCF, 'r') as ifile:
            self.fields = [x for x in ifile.header.info.keys() if x.endswith('_PCTL')]
        self.chrom = None
        self.start = None
        self.stop = None
        self.data = None
    def descriptions(self):
        with pysam.VariantFile(self.VCF, 'r') as ifile:
            for key, value in ifile.header.info.iteritems():
                if key.endswith('_PCTL'):
                     yield '##INFO=<ID={},Number={},Type={},Description="{}">'.format(key, value.number, value.type, value.description)
    def overlaps(self, chrom, position):
        if self.chrom is None or self.start is None or self.stop is None or self.data is None or self.chrom != chrom or self.start > position or self.stop < position:
            return False
        return True
    def get(self, chrom, position, ref, alts):
        if self.has_chr_prefix and not chrom.startswith('chr'):
            chrom = 'chr' + chrom
        elif not self.has_chr_prefix and chrom.startswith('chr'):
            chrom = chrom[3:]
        if not self.overlaps(chrom, position):
            self.chrom = chrom
            self.start = position
            self.stop = position + self.step_bp
            self.data = dict()
            with pysam.VariantFile(self.VCF, 'r') as ifile:
                for record in ifile.fetch(self.chrom, self.start - 1, self.stop + 1):
                    name = ':'.join([record.chrom, str(record.pos), record.ref, ','.join(record.alts)])
                    self.data[name] = [(x, record.info[x]) for x in self.fields]
        return self.data.get(':'.join([chrom, str(position), ref, ','.join(alts)]), [])


def add_percentiles(in_VCF, in_pctl_files, out_VCF):
    percentiles = [Percentiles(x) for x in in_pctl_files]
    with pysam.VariantFile(in_VCF, 'r') as ifile, pysam.BGZFile(out_VCF, 'w') as ofile:
        for p in percentiles:
            for desc in p.descriptions():
                ifile.header.add_line(desc)
        ofile.write('{}'.format(ifile.header))
        for record in ifile:
            for p in percentiles:
                for key, values in p.get(record.chrom, record.pos, record.ref, record.alts):
                    record.info[key] = values
            ofile.write('{}'.format(record))


if __name__ == "__main__":
    args = argparser.parse_args()
    add_percentiles(args.in_VCF, args.in_pctl_files, args.out_VCF)

