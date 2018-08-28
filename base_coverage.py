import json
import time

import pysam
from lookups import IntervalSet
from utils import Xpos


class CoverageHandler(object):
    '''contains coverage (at multiple binning levels) for all chroms'''
    def __init__(self, coverage_files):
        self._single_chrom_coverage_handlers = {}
        for cf in coverage_files:
            coverage_file = CoverageFile(cf['path'], cf.get('binned',False))
            for chrom in coverage_file.get_chroms():
                if chrom not in self._single_chrom_coverage_handlers:
                    self._single_chrom_coverage_handlers[chrom] = SingleChromCoverageHandler(chrom)
                self._single_chrom_coverage_handlers[chrom].add_coverage_file(coverage_file, cf.get('bp-min-length',0))
    def get_coverage_for_intervalset(self, intervalset):
        st = time.time()
        try: single_chrom_coverage_handler = self._single_chrom_coverage_handlers[intervalset.chrom]
        except KeyError: print 'Warning: No coverage for chrom', intervalset.chrom; return []
        coverage = []
        intervalset_length = intervalset.get_length()
        for pair in intervalset.to_obj()['list_of_pairs']:
            coverage.extend(single_chrom_coverage_handler.get_coverage_for_range(pair[0], pair[1], length=intervalset_length))
        print '## COVERAGE: spent {:.3f} seconds tabixing {} coverage bins'.format(time.time()-st, len(coverage))
        return coverage

class SingleChromCoverageHandler(object):
    '''contains coverage (at multiple binning levels) for one chrom'''
    def __init__(self, chrom):
        self._chrom = chrom
        self._coverage_files = []
    def add_coverage_file(self, coverage_file, min_length_in_bases):
        self._coverage_files.append({'coverage_file':coverage_file, 'bp-min-length':min_length_in_bases})
        self._coverage_files.sort(key=lambda d:d['bp-min-length'])
    def get_coverage_for_range(self, start, stop, length=None):
        if length is None: length = stop - start
        assert len(self._coverage_files) >= 1, (self._chrom, start, stop, length, self._coverage_files, str(self))
        assert self._coverage_files[0]['bp-min-length'] <= length, (self._chrom, start, stop, length, self._coverage_files, str(self))
        # get the last (ie, longest `bp-min-length`) coverage_file that has a `bp-min-length` <= length
        coverage_file = next(cf['coverage_file'] for cf in reversed(self._coverage_files) if cf['bp-min-length'] <= length)
        return coverage_file.get_coverage(self._chrom, start, stop)
    def __str__(self):
        return '<SingleChromCoverageHandler chrom={} coverage_files={!r}>'.format(self._chrom, self._coverage_files)
    __repr__ = __str__

class CoverageFile(object):
    '''handles a single tabixed coverage file with any number of chroms'''
    # our coverage files don't include `chr`, so this class prepends `chr` to output and strips `chr` from input
    def __init__(self, path, binned):
        self._tabixfile = pysam.TabixFile(path)
        self._binned = binned
    def get_chroms(self):
        return self._tabixfile.contigs
    def get_coverage(self, chrom, start, stop):
        if not self._binned:
            for row in self._tabixfile.fetch(chrom, start, stop+1, parser=pysam.asTuple()):
                yield json.loads(row[2])
        else:
            # Right now we don't include the region_end column in our coverage files,
            # so there's no way to make sure we get the bin overlapping the start of our query region.
            # To deal with it for now, we'll just use start-50
            # TODO: include region_end in coverage files.
            for row in self._tabixfile.fetch(chrom, max(1, start-50), stop+1, parser=pysam.asTuple()):
                d = json.loads(row[2])
                if d['end'] < start or d['start'] > stop: continue
                d['start'] = max(d['start'], start)
                d['end'] = min(d['end'], stop)
                yield d
    def __str__(self):
        return '<CoverageFile chroms={} path={}>'.format(','.join(self.get_chroms()), self._tabixfile.filename)
    __repr__ = __str__
