import sys
import pysam
import numpy
import re
import gzip
from contextlib import closing
import argparse

argparser = argparse.ArgumentParser(description = 'Computes depth and genotype quality histograms, allele counts, homozygous counts, and etc. Will not remove monomorphic variants.')
argparser.add_argument('--in-gzvcf', metavar = 'file', dest = 'inGZVCF', required = True, help = 'Input VCF file indexed with tabix')
argparser.add_argument('--in-samples', metavar = 'file', dest = 'inSamplesList', default = None, help = 'List of samples to use (one sample per line)')
argparser.add_argument('--chrom', metavar = 'name', dest = 'contig', required = True, help = 'Chromosome name')
argparser.add_argument('--start', metavar = 'bp', dest = 'start', required = True, type=long, help = 'Start position in base-pairs (bp)')
argparser.add_argument('--end', metavar = 'bp', dest = 'end', required = True, type=long, help = 'End position in base-pairs (bp)')
argparser.add_argument('--out-gzvcf', metavar = 'file', dest = 'outGZVCF', required = True, help = 'Output VCF file compressed with gzip')

# Every coverage that is >100 we will put into the last histgram bin (95, 100].
histogram_bins = [0, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51, 56, 61, 66, 71, 76, 81, 86, 91, 96, sys.maxint]
histogram_mids = [2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5, 82.5, 87.5, 92.5, 97.5]

def selectGTFields(inGZVCF, inSamplesList = None):
   inSamples = set()
   inGTFields = []

   if inSamplesList:
      with open(inSamplesList) as f:
         for line in f:
            inSamples.add(line.rstrip())

   with gzip.GzipFile(inGZVCF) as z:
      for line in z:
         if line.startswith('#CHROM'):
            break
      
      if not line:
         raise Error('VCF header was not found!')

      fields = line.rstrip().split('\t')
      
      if inSamplesList:
         for i in xrange(9, len(fields)):
            if fields[i] in inSamples:
               inGTFields.append(i)
      else:
         return range(9, len(fields))
     
   return inGTFields

def processVCF(inGZVCF, inGTFields, contig, start, end, outGZVCF):
   with closing(pysam.Tabixfile(inGZVCF)) as tabix, gzip.GzipFile(outGZVCF, 'w') as out:
      GQIndex = None
      DPIndex = None
      ADIndex = None
      ODIndex = None
      PLIndex = None

      out.write('##INFO=<ID=AVGDP,Number=1,Type=Float,Description="Average Depth per Sample">\n')
      out.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth at Site">\n')
      out.write('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Coverage">\n')
      out.write('##INFO=<ID=AN,Number=1,Type=Integer,Description="Number of Alleles in Samples with Coverage">\n')
      out.write('##INFO=<ID=AC,Number=A,Type=Integer,Description="Alternate Allele Counts in Samples with Coverage">\n')
      out.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Alternate Allele Frequencies">\n')
      out.write('##INFO=<ID=DP_HIST,Number=A,Type=String,Description="Histogram for DP; Mids: %.1f' % histogram_mids[0])
      for i in xrange(1, len(histogram_mids)):
         out.write('|%.1f' % histogram_mids[i])
      out.write('">\n')
      out.write('##INFO=<ID=GQ_HIST,Number=A,Type=String,Description="Histogram for GQ; Mids: %.1f' % histogram_mids[0])
      for i in xrange(1, len(histogram_mids)):
         out.write('|%.1f' % histogram_mids[i])
      out.write('">\n')
      out.write('##INFO=<ID=Hom,Number=A,Type=Integer,Description="Homozygous Counts">\n')
      out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

      for row in tabix.fetch(contig, start - 1, end, parser = pysam.asTuple()):
         altAlleles = map(str, range(1, len(row[4].split(',')) + 1))
         altGQs = dict()
         altDPs = dict()
         altHoms = dict()
         altCounts = dict()
         for altAllele in altAlleles:
            altGQs[altAllele] = []
            altDPs[altAllele] = []
            altHoms[altAllele] = 0
            altCounts[altAllele] = 0      

         fields = row[8].split(':') 
         GQIndex = fields.index('GQ')
         try:
            DPIndex = fields.index('DP')
            ADIndex = None
            ODIndex = None
         except ValueError:
            DPIndex = None
            ADIndex = fields.index('AD')
            ODIndex = fields.index('OD')
         PLIndex = fields.index('PL')

         GQs = []
         DPs = []
         GQ = None
         DP = None
         PL = None
         totalDP = 0
         NS = 0
         AN = 0

         for i in inGTFields:
            fields = row[i].split(':')
            PL = numpy.sum(map(int, fields[PLIndex].split(',')))
            if PL == 0:
               continue

            NS += 1
            AN += 2
        
            GQ = float(fields[GQIndex])
            if DPIndex is None:
               DP = numpy.sum(map(int, fields[ADIndex].split(',')))
               DP += int(fields[ODIndex])
            else:
               DP = int(fields[DPIndex])
         
            GQs.append(GQ)
            DPs.append(DP)
         
            totalDP += DP        

            genoAltAlleles = re.split('/|', fields[0])
            for altAllele in altAlleles:
               count = genoAltAlleles.count(altAllele)
               if count > 0:
                  altGQs[altAllele].append(GQ)
                  altDPs[altAllele].append(DP)
                  if count == 2:
                     altHoms[altAllele] += 1
                  altCounts[altAllele] += count

         GQsHist = numpy.histogram(GQs, bins = histogram_bins)[0]
         DPsHist = numpy.histogram(DPs, bins = histogram_bins)[0]
         altGQsHists = dict()
         altDPsHists = dict()
         for altAllele in altAlleles:
            altGQsHists[altAllele] = numpy.histogram(altGQs[altAllele], bins = histogram_bins)[0]
            altDPsHists[altAllele] = numpy.histogram(altDPs[altAllele], bins = histogram_bins)[0]

         out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t' % (row[0], row[1], row[2], row[3], row[4], row[5], row[6]))
         if DPIndex is None:
            out.write('AVGDP=%.6f' % ((totalDP / float(NS)) if NS > 0 else 0.0)) 
         else:
            out.write('DP=%d;NS=%d' % (totalDP, NS))
         out.write(';AN=%d;AC=%d' % (AN, altCounts[altAlleles[0]]))
         for i in xrange(1, len(altAlleles)):
            out.write(',%d' % altCounts[altAlleles[i]])
         out.write(';AF=%.6f' % ((altCounts[altAlleles[0]] / float(AN)) if AN > 0 else 0.0))
         for i in xrange(1, len(altAlleles)):
            out.write(',%.6f' % ((altCounts[altAlleles[i]] / float(AN)) if AN > 0 else 0.0))
         out.write(';DP_HIST=%d' % DPsHist[0])
         for i in xrange(1, len(DPsHist)):
            out.write('|%d' % DPsHist[i])
         for altAllele in altAlleles:
            out.write(',%d' % altDPsHists[altAllele][0])
            for i in xrange(1, len(altDPsHists[altAllele])):
               out.write('|%d' % altDPsHists[altAllele][i])

         out.write(';GQ_HIST=%d' % GQsHist[0])
         for i in xrange(1, len(GQsHist)):
            out.write('|%d' % GQsHist[i])
         for altAllele in altAlleles:
            out.write(',%d' % altGQsHists[altAllele][0])
            for i in xrange(1, len(altGQsHists[altAllele])):
               out.write('|%d' % altGQsHists[altAllele][i])

         out.write(';Hom=%d' % altHoms['1'])
         for i in xrange(1, len(altAlleles)):
            out.write(',%d' % altHoms[altAlleles[i]])
 
         out.write('\n')

if __name__ == "__main__":
   args = argparser.parse_args()

   inGTFields = selectGTFields(args.inGZVCF, args.inSamplesList)
   processVCF(args.inGZVCF, inGTFields, args.contig, args.start, args.end, args.outGZVCF)

