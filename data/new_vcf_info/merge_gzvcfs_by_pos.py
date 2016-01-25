import sys
import gzip
from collections import deque
import argparse
import operator

argparser = argparse.ArgumentParser(description = 'Merges overlapping VCF (compressed with gzip) files. Requirements: (a) All VCF files must store same chromosome; (b) No regions can be covered by more than two VCF files; (c) There are no such two VCF files, that store nested regions.')
argparser.add_argument('--in', metavar = 'file', dest = 'inGZVCFList', required = True, help = 'List of VCF (compressed with gzip) files. One file per line.')
argparser.add_argument('--out', metavar = 'file', dest = 'outGZVCF', required = True, help = 'Output VCF (compressed with gzip) file')

def readGZVCFsList(inGZVCFsList):
   unorderedGZVCFs = dict()
   with open(inGZVCFsList, 'r') as f:
      for line in f:
         if line.startswith('#'):
            continue
         GZVCF = line.rstrip()
         with gzip.GzipFile(GZVCF) as z:
            for line in z:
               if not line.startswith('#'):
                  fields = line.rstrip().split('\t', 2)
                  unorderedGZVCFs[long(fields[1])] =(GZVCF, fields[0])
                  break

   GZVCFs = [[unorderedGZVCFs[position][0], unorderedGZVCFs[position][1], position, sys.maxsize] for position in sorted(unorderedGZVCFs.iterkeys())]   
 
   for i in xrange(1, len(GZVCFs)):
      if GZVCFs[i - 1][1] != GZVCFs[i][1]:
         raise Exception('Input files store different contigs/chromosomes!')
      if GZVCFs[i - 1][2] == GZVCFs[i][2]:
         raise Exception('Two input files store identical leftmost positions!')

   for i in xrange(0, len(GZVCFs) - 1):
      GZVCFs[i][3] = GZVCFs[i + 1][2] 

   return GZVCFs

def mergeGZVCFsByPos(GZVCFs, outGZVCF):
   with gzip.GzipFile(outGZVCF, 'w') as out:

      with gzip.GzipFile(GZVCFs[0][0]) as z:
         for line in z:
            if line.startswith('#'):
               out.write(line)
            else:
               break

      overlapHead = deque([])
      overlapTail = deque([])
  
      for GZVCF in GZVCFs:
         with gzip.GzipFile(GZVCF[0]) as z: 
            for line in z:
               if line.startswith('#'):
                  continue
               line = line.rstrip()
               if not line:
                  continue

               fields = line.split('\t', 2)

               contig = fields[0]
               if contig != GZVCF[1]:
                  raise Exception('Multiple chromosomes detected within single VCF!')

               position = long(fields[1])

               while overlapHead:
                  (overlapPosition, overlapLine) = overlapHead[0]

                  if overlapPosition >= GZVCF[3]:
                     raise Exception("Overlapping regions are present in more than two VCF files!")

                  if overlapPosition < position:
                     out.write(overlapLine + '\n')
                     overlapHead.popleft()
                  elif overlapPosition == position:
                     if overlapLine != line:
                        out.write(overlapLine + '\n')
                     overlapHead.popleft()
                  else:
                     break

               if position < GZVCF[3]:
                  out.write(line + '\n')
               else:
                  overlapTail.append((position, line))

         if overlapHead:
            raise Exception("Nested regions detected in two VCF files!")
      
         overlapHead = overlapTail
         overlapTail = deque([])
     
      if overlapTail:
         raise Exception("Error while merging VCF files!")  


if __name__ == "__main__":
   args = argparser.parse_args()

   GZVCFs = readGZVCFsList(args.inGZVCFList)
   mergeGZVCFsByPos(GZVCFs, args.outGZVCF)

