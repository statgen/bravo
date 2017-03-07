import gzip
import json
import sys
import pysam
import argparse
import numpy
import collections
from contextlib import closing

argparser = argparse.ArgumentParser(description = 'Aggregate depth information (output as JSON) from individual depth files (generated using SAMtools mpileup).')
argparser.add_argument('--in', metavar = 'file', dest = 'inFileList', required = True, help = 'Input file listing all depth files (one depth file per sample) generated using SAMtools mpileup. One file per line.')
argparser.add_argument('--chromosome', metavar = 'name', dest = 'chromosome', required = True, help = 'Chromosome name.')
argparser.add_argument('--start', metavar = 'bp', dest = 'startbp', type = long, required = True, help = 'Region start position in bp.')
argparser.add_argument('--end', metavar = 'bp', dest = 'endbp', type = long, required = True, help = 'Region end position in bp.')

breaks = [1, 5, 10, 15, 20, 25, 30, 50, 100]

def readDepthChunk(depthFile, contig, start, end):
   chunk = dict()
   with closing(pysam.Tabixfile(depthFile)) as tabix:
      for row in tabix.fetch(contig, start - 1, end, parser = pysam.asTuple()):
         chunk[long(row[1])] = int(row[3])
   return chunk

def readDepthFileList(inFileList):
   depthFiles = list()
   with open(inFileList, 'r') as ifile:
      for line in ifile:
         line = line.rstrip()
         if not line:
            continue
         depthFiles.append(line)
   return depthFiles

def readDepthChunks(depthFiles, contig, start, end):
   chunks = dict()
   for depthFile in depthFiles:
      chunk = readDepthChunk(depthFile, contig, start, end)
      for position, depth in chunk.iteritems():
         if position in chunks:
            chunks[position].append(depth)
         else:
            chunks[position] = [depth]
   return chunks

def writeDepthChunks(chunks, chromosome, nDepthFiles):
   for position, depths in collections.OrderedDict(sorted(chunks.items())).iteritems():
      counts = [0] * len(breaks) 
      for depth in depths:
         for i in xrange(len(breaks) - 1, -1, -1):
            if depth >= breaks[i]:
               counts[i] += 1
               break
      for i in xrange(len(breaks) - 2, -1, -1):
         counts[i] += counts[i + 1];

      sys.stdout.write('%s\t%d\t{"chrom":"%s","start":%d,"end":%d,"mean":%g,"median":%g' % (chromosome, position, chromosome, position, position, numpy.mean(depths), numpy.median(depths)))
      for i in xrange(0, len(breaks)): 
         sys.stdout.write(',"%d":%g' % (breaks[i], counts[i] / nDepthFiles))
      sys.stdout.write('}\n')

if __name__ == '__main__':
   args = argparser.parse_args()
   depthFiles = readDepthFileList(args.inFileList)
   writeDepthChunks(readDepthChunks(depthFiles, args.chromosome, args.startbp, args.endbp), args.chromosome, float(len(depthFiles)))

