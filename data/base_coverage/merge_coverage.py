import sys
import pysam
import gzip
import json
from collections import deque
import argparse
import operator


argparser = argparse.ArgumentParser(description = 'Merges overlapping JSON coverage (compressed with bgzip/gzip) files. Requirements: (a) All files must store same chromosome; (b) No regions can be covered by more than two files; (c) There are no such two files, that store nested regions; (d) All positions within single coverage file are unique and in ascending order.')
argparser.add_argument('-i', '--in', metavar = 'file', dest = 'in_files_list', required = True, help = 'List of JSON coverage (compressed with bgzip/gzip) files. One file per line.')
argparser.add_argument('-o', '--out', metavar = 'file', dest = 'out_merged_file', required = True, help = 'Output JSON coverage (compressed with bgzip) file')


def read_files_list(files_list):
   unordered_coverage_files = dict()
   with open(files_list, 'r') as ifile:
      for line in ifile:
         if line.startswith('#'):
            continue
         coverage_file = line.rstrip()
         with gzip.GzipFile(coverage_file) as iz:
            for line in iz:
               fields = line.rstrip().split('\t', 2)
               position = long(fields[1])
               chrom = fields[0]
               unordered_coverage_files[position] = (coverage_file, chrom)
               break
   ordered_coverage_files = [{
            'name': unordered_coverage_files[position][0],
            'chrom': unordered_coverage_files[position][1],
            'leftmost_position': position,
            'next_leftmost_position': sys.maxsize} for position in sorted(unordered_coverage_files.iterkeys())]
   for i in xrange(1, len(ordered_coverage_files)):
      if ordered_coverage_files[i - 1]['chrom'] != ordered_coverage_files[i]['chrom']:
         raise Exception('Input files store different contigs/chromosomes!')
      if ordered_coverage_files[i - 1]['leftmost_position'] == ordered_coverage_files[i]['leftmost_position']:
         raise Exception('Two input files store identical leftmost positions!')
   for i in xrange(0, len(ordered_coverage_files) - 1):
      ordered_coverage_files[i]['next_leftmost_position'] = ordered_coverage_files[i + 1]['leftmost_position']
   return ordered_coverage_files


def merge_coverage_files(coverage_files, out_coverage_file):
   with pysam.BGZFile(out_coverage_file, 'w') as oz:
      overlap_head = deque([])
      overlap_tail = deque([])
      for coverage_file in coverage_files:
         last_position = None
         with gzip.GzipFile(coverage_file['name']) as iz:
            for line in iz:
               fields = line.split('\t')
               chrom = fields[0]
               if chrom != coverage_file['chrom']:
                  raise Exception('Multiple chromosomes detected within {} coverage file!'.format(coverage_file['name']))
               position = long(fields[1])

               if last_position is None or last_position < position:
                  last_position = position
               else:
                  raise Exception('Positions within {} coverage file are not in ascending order or not unique!'.format(coverage_file['name']))

               while overlap_head:
                  (overlap_position, overlap_line) = overlap_head[0]

                  if overlap_position >= coverage_file['next_leftmost_position']:
                     raise Exception("Overlapping regions are present in more than two coverage files!")

                  if overlap_position < position:
                     oz.write(overlap_line)
                     overlap_head.popleft()
                  elif overlap_position == position:
                     overlap_head.popleft()
                     overlap_data = json.loads(overlap_line.split('\t')[2])
                     data = json.loads(fields[2])
                     if overlap_data['mean'] > data['mean']:
                        line = overlap_line
                     else:
                        break
                  else:
                     break

               if position < coverage_file['next_leftmost_position']:
                  oz.write(line)
               else:
                  overlap_tail.append((position, line))

         if overlap_head:
            raise Exception("Nested regions detected in two coverage files!")
         overlap_head = overlap_tail
         overlap_tail = deque([])

      if overlap_tail:
         raise Exception("Error while merging coverage files!")


if __name__ == "__main__":
   args = argparser.parse_args()
   coverage_files = read_files_list(args.in_files_list)
   merge_coverage_files(coverage_files, args.out_merged_file)
   pysam.tabix_index(args.out_merged_file, seq_col = 0, start_col = 1, end_col = 1)
