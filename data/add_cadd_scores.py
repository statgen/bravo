import sys
import gzip
import pysam
import re
from contextlib import closing
import argparse

argparser = argparse.ArgumentParser(description = 'Adds CADD scores to VCF.')
argparser.add_argument('-i', '--in-vcf', metavar = 'file', dest = 'in_VCF', required = True, help = 'Input VCF/BCF file. File may be compressed with gzip/bgzip.')
argparser.add_argument('-c', '--in-cadd', metavar = 'file', dest = 'in_CADD_files', nargs = '+', required = True, help = 'Input CADD (from http://cadd.gs.washington.edu) files. Must be compressed with bgzip and indexed with tabix.')
argparser.add_argument('-o', '--out-vcf', metavar = 'file', dest = 'out_VCF', required = True, help = 'Output VCF file. File will be compressed with bgzip.')


class CADD(object):
    def __init__(self, files):
        self.step_bp = 1000000
        self.files = files
        self.has_chr_prefix = False
        for f in self.files:
            with pysam.Tabixfile(f, 'r') as tabix:
                if any(c.startswith('chr') for c in tabix.contigs):
                    self.has_chr_prefix = True
                    break
        self.chrom = None
        self.start = None
        self.stop = None
        self.data = None
    def overlaps(self, chrom, position):
        if self.chrom is None or self.start is None or self.stop is None or self.data is None or self.chrom != chrom or self.start > position or self.stop < position:
            return False
        return True
    def get(self, chrom, position, ref, alt):
        if self.has_chr_prefix and not chrom.startswith('chr'):
            chrom = 'chr' + chrom
        elif not self.has_chr_prefix and chrom.startswith('chr'):
            chrom = chrom[3:]
        if not self.overlaps(chrom, position):
            self.chrom = chrom
            self.start = position
            self.stop = position + self.step_bp
            self.data = dict()
            for f in self.files:
                with pysam.Tabixfile(f, 'r') as tabix:
                    for row in tabix.fetch(self.chrom, self.start - 1, self.stop + 1, parser = pysam.asTuple()):
                        name = ':'.join(row[:4])
                        cadd_raw, cadd_phred = map(float, row[4:6])
                        if name in self.data:
                            if self.data[name][1] < cadd_phred:
                                self.data[name] = (cadd_raw, cadd_phred)
                        else:
                            self.data[name] = (cadd_raw, cadd_phred)
        return self.data.get(':'.join((chrom, str(position), ref, alt)), (None, None))


def addCADD(in_VCF, in_CADD_files, out_VCF):
    cadds = CADD(in_CADD_files)
    with pysam.VariantFile(in_VCF, 'r') as ifile, pysam.BGZFile(out_VCF, 'w') as ofile:
        for x in ['CADD_RAW', 'CADD_PHRED']:
            if x in ifile.header.info:
                raise Exception('{} already exists in input VCF/BCF.'.format(x))
        ifile.header.add_line('##INFO=<ID=CADD_RAW,Number=A,Type=Float,Description="Raw CADD scores">')
        ifile.header.add_line('##INFO=<ID=CADD_PHRED,Number=A,Type=Float,Desctiption="Phred-scaled CADD scores">')
        ofile.write('{}'.format(ifile.header))
        for record in ifile:
            raw_scores = []
            phred_scores = []
            for alt in record.alts:
                cadd_raw, cadd_phred = cadds.get(record.chrom, record.pos, record.ref, alt)
                raw_scores.append(cadd_raw)
                phred_scores.append(cadd_phred)
            if any(x for x in raw_scores) and any(x for x in phred_scores):
                record.info['CADD_RAW'] = raw_scores
                record.info['CADD_PHRED'] = phred_scores
            ofile.write('{}'.format(record))


if __name__ == "__main__":
   args = argparser.parse_args()
   addCADD(args.in_VCF, args.in_CADD_files, args.out_VCF)

