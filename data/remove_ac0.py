import sys
import gzip
import argparse

argparser = argparse.ArgumentParser(description = 'Removes variants that have at least one allele count (AC) equal to 0 or AN.')
argparser.add_argument('--in', metavar = 'file', dest = 'inGZVCF', required = True, help = 'Input VCF file, compressed with gzip/bgzip.')
argparser.add_argument('--out', metavar = 'file', dest = 'outGZVCF', required = True, help = 'Ouput VCF file, compressed with gzip.')

def removeAC0(inGZVCF, outGZVCF):
   with gzip.GzipFile(inGZVCF, 'r') as iz, gzip.GzipFile(outGZVCF, 'w') as oz:
      for line in iz:
         if line.startswith('#'):
            oz.write(line)
            continue

         info = line.rstrip().split('\t')[7]
         ac = None
         an = None

         for field in info.split(';'):
            key, value = field.split('=', 1)
            if key == 'AC':
               ac = map(int, value.split(','))
            elif key == 'AN':
               an = int(value)

         if not ac:
            raise Exception('AC INFO fields was not found!')

         if not an:
            raise Exception('AN INFO fields was not found!')

         if any(an == value or value == 0 for value in ac):
            continue

         oz.write(line)
 
if __name__ == '__main__':
   args = argparser.parse_args()
   removeAC0(args.inGZVCF, args.outGZVCF)
