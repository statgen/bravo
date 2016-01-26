import sys
import gzip
import json
from collections import deque
import argparse
import operator

argparser = argparse.ArgumentParser(description = 'Merges overlapping coverage (compressed with gzip) files. Requirements: (a) All files must store same chromosome; (b) No regions can be covered by more than two files; (c) There are no such two files, that store nested regions; (d) All positions within single coverage file are unique and in ascending order.')
argparser.add_argument('--in', metavar = 'file', dest = 'inCoverageGZFsList', required = True, help = 'List of coverage (compressed with gzip) files. One file per line.')
argparser.add_argument('--out', metavar = 'file', dest = 'outCoverageGZF', required = True, help = 'Output coverage (compressed with gzip) file')

def readCoverageGZFsList(inCoverageGZFsList):
   unorderedCoverageGZFs = dict()
   with open(inCoverageGZFsList, 'r') as f:
      for line in f:
         if line.startswith('#'):
            continue
         coverageGZF = line.rstrip()
         with gzip.GzipFile(coverageGZF) as z:
            for line in z:
               fields = line.rstrip().split('\t', 2)
               unorderedCoverageGZFs[long(fields[1])] =(coverageGZF, fields[0])
               break

   coverageGZFs = [[unorderedCoverageGZFs[position][0], unorderedCoverageGZFs[position][1], position, sys.maxsize] for position in sorted(unorderedCoverageGZFs.iterkeys())]   
 
   for i in xrange(1, len(coverageGZFs)):
      if coverageGZFs[i - 1][1] != coverageGZFs[i][1]:
         raise Exception('Input files store different contigs/chromosomes!')
      if coverageGZFs[i - 1][2] == coverageGZFs[i][2]:
         raise Exception('Two input files store identical leftmost positions!')

   for i in xrange(0, len(coverageGZFs) - 1):
      coverageGZFs[i][3] = coverageGZFs[i + 1][2] 

   return coverageGZFs

def mergeCoverageGZFs(coverageGZFs, outCoverageGZF):
   with gzip.GzipFile(outCoverageGZF, 'w') as oz:

      overlapHead = deque([])
      overlapTail = deque([])
 
      for coverageGZF in coverageGZFs:
         lastPosition = None
         with gzip.GzipFile(coverageGZF[0]) as iz: 
            for line in iz:
               fields = line.split('\t')

               contig = fields[0]
               if contig != coverageGZF[1]:
                  raise Exception('Multiple chromosomes detected within single coverage file!')
               position = long(fields[1])
              
               if lastPosition is None or lastPosition < position:
                  lastPosition = position
               else:
                  raise Exception('Positions within single coverage file are not in ascending order or not unique!')
 
               while overlapHead:
                  (overlapPosition, overlapLine) = overlapHead[0]

                  if overlapPosition >= coverageGZF[3]:
                     raise Exception("Overlapping regions are present in more than two coverage files!")

                  if overlapPosition < position:
                     oz.write(overlapLine)
                     overlapHead.popleft()
                  elif overlapPosition == position:
                     overlapHead.popleft()
                     overlapData = json.loads(overlapLine.split('\t')[2])
                     data = json.loads(fields[2])
                     if overlapData['mean'] > data['mean']:
                        line = overlapLine
                     else:
                        break 
                  else:
                     break

               if position < coverageGZF[3]:
                  oz.write(line)
               else:
                  overlapTail.append((position, line))
                 
         if overlapHead:
            raise Exception("Nested regions detected in two coverage files!")
      
         overlapHead = overlapTail
         overlapTail = deque([])
     
      if overlapTail:
         raise Exception("Error while merging coverage files!")  


if __name__ == "__main__":
   args = argparser.parse_args()

   coverageGZFs = readCoverageGZFsList(args.inCoverageGZFsList)
   mergeCoverageGZFs(coverageGZFs, args.outCoverageGZF)

