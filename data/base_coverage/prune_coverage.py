import gzip
import json
import sys
import argparse

argparser = argparse.ArgumentParser(description = 'Prunes base coverage by grouping bases into bins of similar median and mean depths. Starting with base X that has mean and median depths Z(Y) and Y(X), every next base X+1 is added to the same bin if |Z(X+1) - Z(X)| <= LIMIT and |Y(X+1) - Y(X)| <= LIMIT')
argparser.add_argument('--in', metavar = 'file', dest = 'inCoverageGZF', required = True, help = 'Input coverage (compressed with gzip/bgzip) file')
argparser.add_argument('--limit', metavar = 'float', dest = 'fluctuationLimit', required = True, type = float, help = 'Threshold for maximal fluctuation of median and mean depths within a bin.')
argparser.add_argument('--out', metavar = 'file', dest = 'outCoverageGZF', required = True, help = 'Output coverage (compressed with gzip) file')

def writeData(outputFile, jsonData):
   outputFile.write(jsonData['chrom'])
   outputFile.write('\t')
   outputFile.write(str(jsonData['start']))
   outputFile.write('\t')
   json.dump(jsonData, outputFile, separators = (',', ':'))
   outputFile.write('\n')            

def prune(inCoverageGZF, outCoverageGZF, fluctuationLimit = 0.25):
   with gzip.GzipFile(inCoverageGZF, 'r') as iz, gzip.GzipFile(outCoverageGZF, 'w') as oz:
      fields = iz.readline().split('\t')
      bin_data = json.loads(fields[2])
      for line in iz:
         fields = line.split('\t')
         data = json.loads(fields[2])
         if ((abs(bin_data['mean'] - data['mean']) > fluctuationLimit) or (abs(bin_data['median'] - data['median']) > fluctuationLimit)):
             writeData(oz, bin_data)
             bin_data = data
         bin_data['end'] = data['start']
      writeData(oz, bin_data)

if __name__ == "__main__":
   args = argparser.parse_args()
   
   prune(args.inCoverageGZF, args.outCoverageGZF, args.fluctuationLimit)
