import gzip
import os
import json
import sys
import pysam
import argparse
import numpy
import collections
from contextlib import closing

argparser = argparse.ArgumentParser(description = 'Aggregate depth information (output as JSON) from individual depth files (generated using SAMtools mpileup).')
argparser.add_argument('-i', '--in', metavar = 'file', dest = 'inFileList', required = True, help = 'Input file listing all depth files (one depth file per sample) generated using SAMtools mpileup. One file per line.')

sub_argparsers = argparser.add_subparsers(dest = 'subparser_name')

chunk_argparser = sub_argparsers.add_parser('chunk', help = 'Generate and print chromosome chunks.')
chunk_argparser.add_argument('-c', '--chromosome', metavar = 'name', dest = 'chromosome', required = True, help = 'Chromosome name.')
chunk_argparser.add_argument('-s', '--size', metavar = 'bp', dest = 'chunk_size_bp', type = long, required = True, help = 'Maximal chunk size in base-pairs.')

aggregate_argparser = sub_argparsers.add_parser('aggregate', help = 'Aggregate depth information across individuals.')
aggregate_argparser.add_argument('-c', '--chromosome', metavar = 'name', dest = 'chromosome', required = True, help = 'Chromosome name.')
aggregate_argparser.add_argument('-s', '--start', metavar = 'bp', dest = 'startbp', type = long, required = True, help = 'Region start position in bp.')
aggregate_argparser.add_argument('-e', '--end', metavar = 'bp', dest = 'endbp', type = long, required = True, help = 'Region end position in bp.')

breaks = [1, 5, 10, 15, 20, 25, 30, 50, 100]

def getMinMaxPositions(depthFile, contig):
    with closing(pysam.TabixFile(depthFile)) as tabix:
        first_entry = None
        for first_entry in tabix.fetch(contig, 0, parser = pysam.asTuple()):
            break
        last_Mbp = 0
        while any(True for _ in tabix.fetch(contig, last_Mbp, parser = pysam.asTuple())):
            last_Mbp += 5000000
        last_entry = None
        for last_entry in tabix.fetch(contig, last_Mbp - 5000000, parser = pysam.asTuple()):
            pass
        return (long(first_entry[1]) if first_entry is not None else None, long(last_entry[1]) if last_entry is not None else None)

def chunk(inFileList, contig, chunk_size_bp):
    depthFiles = readDepthFileList(args.inFileList)
    starts = []
    ends = []
    for depthFile in depthFiles:
        start, end = getMinMaxPositions(depthFile, contig)
        if start is not None:
            starts.append(start)
        if end is not None:
            ends.append(end)
    start = min(starts)
    end = max(ends)
    n_chunks = int(numpy.ceil((end - start) / float(chunk_size_bp)))
    for i, s in enumerate(xrange(start, end, chunk_size_bp)):
        e = end if i == n_chunks - 1 else s + chunk_size_bp - 1
        output_file = contig + '_' + str(s) + '_' + str(e) + '.json.bgz'
        print 'python', os.path.realpath(__file__), '-i', inFileList, 'aggregate', '-c', contig, '-s', s, '-e', e, '| bgzip -c >', output_file

def readDepthChunk(depthFile, contig, start, end):
   chunk = dict()
   if start > 0:
       start -= 1
   with closing(pysam.Tabixfile(depthFile)) as tabix:
      for row in tabix.fetch(contig, start, end, parser = pysam.asTuple()):
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
   if chromosome.startswith('chr'):
       chromosome = chromosome[3:]
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
   if args.subparser_name == 'chunk':
      chunk(args.inFileList, args.chromosome, args.chunk_size_bp)
   elif args.subparser_name == 'aggregate':
      depthFiles = readDepthFileList(args.inFileList)
      writeDepthChunks(readDepthChunks(depthFiles, args.chromosome, args.startbp, args.endbp), args.chromosome, float(len(depthFiles)))
