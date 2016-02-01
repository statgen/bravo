import sys
import gzip
import re
import argparse

argparser = argparse.ArgumentParser(description = 'Export INFO fields from one VCF file into another VCF file. Both VCF files must store same number of identical variants (i.e. same chromosomes, same genetic positions).')
argparser.add_argument('--import-gzvcf', metavar = 'file', dest = 'importGZVCF', required = True, help = 'Input VCF file compressed with gzip/bgzip. This file will import new INFO fields.')
argparser.add_argument('--export-gzvcf', metavar = 'file', dest = 'exportGZVCF', required = True, help = 'Input VCF file compressed with gzip/bgzip. This file will export its INFO fields')
argparser.add_argument('--out-gzvcf', metavar = 'file', dest = 'outGZVCF', required = True, help = 'Output VCF file compressed with gzip.')
argparser.add_argument('--info', metavar = 'name', dest = 'infoNames', nargs = '+', required = True, help = 'List of INFO fields to import.')
argparser.add_argument('--meta', metavar = 'name', dest = 'metaKeys', nargs = '+', required = False, help = 'List of META keys to import.')

def importINFO(importGZVCF, exportGZVCF, outGZVCF, infoNames, metaKeys):
   metaInfoPattern = re.compile('^##INFO=<ID=([a-zA-Z0-9_\-]+),Number=(A|R|G|\.|[0-9]+),Type=(Integer|Float|Flag|Character|String),Description="((?:[^"]|\\\\")*)"(?:,Source="(?:[^"]|\\\\")*")?(?:,Version="(?:[^"]|\\\\")*")?>')
   metaKeyValuePattern = re.compile('^##([_a-zA-Z0-9]+)=.+$')

   with gzip.GzipFile(importGZVCF, 'r') as iz, gzip.GzipFile(exportGZVCF, 'r') as ez, gzip.GzipFile(outGZVCF, 'w') as oz:
      exportedMetaInfoLines = dict()
      exportedMetaKeyValueLines = dict()
      
      for eline in ez:
         if eline.startswith('##'):
            match = metaKeyValuePattern.match(eline)
            if match:
               if match.group(1) == 'INFO':
                  match = metaInfoPattern.match(eline)
                  if match:
                     if match.group(1) in infoNames:
                        exportedMetaInfoLines[match.group(1)] = eline
                  else:
                     raise Exception('Error while parsing INFO meta-information line in %s!' % exportGZVCF)
               elif metaKeys:
                  if match.group(1) in metaKeys:
                     if match.group(1) in exportedMetaKeyValueLines:
                        exportedMetaKeyValueLines[match.group(1)].append(eline)
                     else:
                        exportedMetaKeyValueLines[match.group(1)] = [eline]
            else:
               raise Exception('Error while parsing meta-information line in %s!' % exportGZVCF)
         else:
            break

      if len(exportedMetaInfoLines) != len(infoNames):
         raise Exception('Some of the specified INFO fields are missing in the %s!' % exportGZVCF)

      if metaKeys and len(exportedMetaKeyValueLines) != len(metaKeys):
         raise Exception('Some of the specified META keys are missing in the %s!' % exportGZVCF)
   
      for iline in iz:
         if iline.startswith('##'):
            match = metaKeyValuePattern.match(iline)
            if match:
               if match.group(1) == 'INFO':
                  match = metaInfoPattern.match(iline)
                  if match:
                     if match.group(1) in infoNames:
                        raise Exception('%s INFO field exists in both VCF files!' % match.group(1))
                  else:
                     raise Exception('Error while parsing INFO meta-information line in %s!' % importGZVCF)
               elif metaKeys:
                  if match.group(1) in exportedMetaKeyValueLines:
                     raise Exception('%s META key exists in both VCF files!' % match.group(1))
            else:
               raise Exception('Error while parsing meta-information line in %s!' % importGZVCF)
            oz.write(iline)
         else:
            break

      for metaKey, lines in exportedMetaKeyValueLines.iteritems():
         for line in lines:
            oz.write(line)

      for metaInfo, line in exportedMetaInfoLines.iteritems():
         oz.write(line)

      oz.write(iline)

      eline = ez.readline()
      iline = iz.readline()
      info = dict()
      while eline and iline:
         efields = eline.rstrip().split('\t')
         ifields = iline.rstrip().split('\t')
         info.clear()       

         if efields[0] != ifields[0] or efields[1] != efields[1]:
            raise Exception('VCF files store different chromosomes and/or genetic positions!')

         for infoField in efields[7].split(';'):
            name, value = infoField.split('=', 1)
            info[name] = value

         oz.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (ifields[0], ifields[1], ifields[2], ifields[3], ifields[4], ifields[5], ifields[6], ifields[7]))
         for infoName in infoNames:
            oz.write(';%s=%s' % (infoName, info[infoName]))
         oz.write('\n')
         eline = ez.readline()
         iline = iz.readline()

      if eline or iline:
         raise Exception('VCF files have different number of variant lines!')

if __name__ == "__main__":
   args = argparser.parse_args()
   importINFO(args.importGZVCF, args.exportGZVCF, args.outGZVCF, args.infoNames, args.metaKeys)
