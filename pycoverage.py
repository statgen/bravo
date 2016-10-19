import json
import pysam
from intervaltree import IntervalTree
import sys
import utils

# TODO:
# - organize first by contig, then by lengths.
# - assert that ranges are (1) never overlapping and (2) cover [0,inf]
#    - so, maybe each (filepath,contig) should just have a min_length.  Then we use bisect.

class Coverage(object):

    def __init__(self):
        self.tabix_path = {}
        self.tabix_handler = {}

    def setTabixPath(self, contig, path):
        self.tabix_path[contig] = path

    def open(self):
        self.tabix_handler = {contig: pysam.Tabixfile(path) for contig, path in self.tabix_path.iteritems()}

    def close(self):
        for contig, handler in self.tabix_handler.iteritems():
            handler.close()

    def getCoverage(self, contig, start, end):
        try:
            handler = self.tabix_handler[contig]
        except KeyError:
            return None
        return [json.loads(row[2]) for row in handler.fetch(contig, start, end + 1, parser = pysam.asTuple())]

    def getCoverageX(self, xstart, xend):
        contig = utils.CHROMOSOME_NUMBER_TO_CODE[xstart // int(1e9)]
        if contig.startswith('chr'): contig = contig[3:] # TODO: this is gross, why do I need to do this?
        start = xstart % int(1e9)
        end = xend % int(1e9)
        return self.getCoverage(contig, start, end)


class CoverageCollection(object):

    def __init__(self):
        self.collection = IntervalTree()

    def setTabixPath(self, min_length, max_length, contig, path):
        coverage = Coverage()
        coverage.setTabixPath(contig, path)
        self.collection[min_length:max_length] = coverage

    def openAll(self):
        for coverage in self.collection.items():
           coverage.data.open()

    def closeAll(self):
        for coverage in self.collection.items():
            coverage.data.close()

    def clearAll(self):
        self.closeAll()
        collection.clear()

    def getCoverage(self, contig, start, end):
        for coverage in self.collection.search(end - start):
            rv = coverage.data.getCoverage(contig, start, end)
            if rv is not None: return rv
        return []

    def getCoverageX(self, xstart, xend):
        for coverage in self.collection.search(xend - xstart):
            rv = coverage.data.getCoverageX(xstart, xend)
            if rv is not None:
                return rv
        return []

    def printAll(self):
        for coverage in self.collection.items():
            print coverage
