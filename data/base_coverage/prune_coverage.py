import gzip
import pysam
import json
import sys
import argparse

argparser = argparse.ArgumentParser(description = 'Prunes base coverage by grouping bases into bins of similar median and mean depths. Starting with base X that has mean and median depths Z(X) and Y(X), every next base X+1 is added to the same bin if |Z(X+1) - Z(X)| <= LIMIT and |Y(X+1) - Y(X)| <= LIMIT')
argparser.add_argument('-i', '--in', metavar = 'file', dest = 'in_coverage_file', required = True, help = 'Input JSON coverage (compressed with gzip/bgzip) file')
argparser.add_argument('-l', '--limit', metavar = 'float', dest = 'fluctuation_limit', required = True, type = float, help = 'Threshold for maximal fluctuation of median and mean depths within a bin.')
argparser.add_argument('-o', '--out', metavar = 'file', dest = 'out_coverage_file', required = True, help = 'Output JSON coverage (compressed with bgzip) file')


def write_data(of, json_data):
   of.write('{}\t{}\t'.format(json_data['chrom'], json_data['start']))
   json.dump(json_data, of, separators = (',', ':'))
   of.write('\n')


def prune(in_coverage_file, out_coverage_file, fluctuation_limit = 0.25):
   with gzip.GzipFile(in_coverage_file, 'r') as iz, pysam.BGZFile(out_coverage_file, 'w') as oz:
      fields = iz.readline().split('\t')
      bin_data = json.loads(fields[2])
      for line in iz:
         fields = line.split('\t')
         data = json.loads(fields[2])
         if ((abs(bin_data['mean'] - data['mean']) > fluctuation_limit) or (abs(bin_data['median'] - data['median']) > fluctuation_limit)):
             write_data(oz, bin_data)
             bin_data = data
         bin_data['end'] = data['start']
      write_data(oz, bin_data)


if __name__ == "__main__":
   args = argparser.parse_args()
   prune(args.in_coverage_file, args.out_coverage_file, args.fluctuation_limit)
   pysam.tabix_index(args.out_coverage_file, seq_col = 0, start_col = 1, end_col = 1)
