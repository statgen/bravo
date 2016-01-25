import sys
import gzip
import pysam
import re
from contextlib import closing
import argparse

argparser = argparse.ArgumentParser(description = 'Adds CADD scores to VCF.')
argparser.add_argument('--in-gzvcf', metavar = 'file', dest = 'inGZVCF', required = True, help = 'Input VCF file compressed with gzip')
argparser.add_argument('--in-cadd', metavar = 'file', dest = 'inCADDFiles', nargs = '+', required = True, help = 'Input CADD (from http://cadd.gs.washington.edu) files indexed with tabix')
argparser.add_argument('--out-gzvcf', metavar = 'file', dest = 'outGZVCF', required = True, help = 'Output VCF file compressed with gzip')

def readCADD(inCADDFiles, contig, start, end):
   CADD = dict()

   for CADDFile in inCADDFiles:
      with closing(pysam.Tabixfile(CADDFile)) as tabix:
         for row in tabix.fetch(contig, start - 1, end, parser = pysam.asTuple()):
            name = row[0] + ':' + row[1] + '_' + row[2] + '/' + row[3]
            cadd_raw = float(row[4])
            cadd_phred = float(row[5])

            if name in CADD:
               if CADD[name][1] < cadd_phred:
                  CADD[name] = (cadd_raw, cadd_phred)
            else:
               CADD[name] = (cadd_raw, cadd_phred)

   return CADD

def addCADD(inGZVCF, inCADDFiles, outGZVCF):
   metaInfoPattern = re.compile('^##INFO=<ID=([a-zA-Z0-9_\-]+),Number=(A|R|G|\.|[0-9]+),Type=(Integer|Float|Flag|Character|String),Description="((?:[^"]|\\\\")*)"(?:,Source="(?:[^"]|\\\\")*")?(?:,Version="(?:[^"]|\\\\")*")?>')

   CADDcontig = ''
   CADDposition = -1
   CADD = None
   stepSize = 1000000

   with gzip.GzipFile(inGZVCF, 'r') as iz, gzip.GzipFile(outGZVCF, 'w') as oz:
      for line in iz:
         if line.startswith('##'):
            match = metaInfoPattern.match(line.rstrip())
            if match: 
               if match.group(1) == 'CADD_RAW':
                  raise Exception('CADD_RAW already exists in input VCF')
               elif match.group(1) == 'CADD_PHRED':
                  raise Exception('CADD_PHRED already exists in input VCF')
            oz.write(line)
         else:
            break
           
      oz.write('##INFO=<ID=CADD_RAW,Number=A,Type=Float,Description="Raw CADD scores">\n')
      oz.write('##INFO=<ID=CADD_PHRED,Number=A,Type=Float,Desctiption="Phred-scaled CADD scores">\n')
      oz.write(line)

      for line in iz:
         if line.startswith('#'):
            oz.write(line)
            continue
         scores = []
         fields = line.rstrip().split('\t')
         contig = fields[0]
         position = long(fields[1])
         if CADDcontig != contig or CADDposition < position:
            CADDcontig = contig
            CADDposition = position + stepSize
            CADD = readCADD(inCADDFiles, CADDcontig, position, CADDposition)
 
         for altAllele in fields[4].split(','):
            name = fields[0] + ':' + fields[1] + '_' + fields[3] + '/' + altAllele 
            if name in CADD:
               scores.append(CADD[name])
         
         oz.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], fields[6], fields[7]))
         if scores:
            oz.write(';CADD_RAW=%.6f' % scores[0][0])
            for i in xrange(1, len(scores)):
               oz.write(',%.6f' % scores[i][0])
            oz.write(';CADD_PHRED=%.6f' % scores[0][1])
            for i in xrange(1, len(scores)):
               oz.write(',%.6f' % scores[i][1])
         oz.write('\n')


if __name__ == "__main__":
   args = argparser.parse_args()
   addCADD(args.inGZVCF, args.inCADDFiles, args.outGZVCF)

