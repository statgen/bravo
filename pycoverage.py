import json
import pysam
from intervaltree import IntervalTree
import sys   
 
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
        return [json.loads(row[2]) for row in self.tabix_handler[contig].fetch(contig, start, end + 1, parser = pysam.asTuple())]

    def getCoverageX(self, xstart, xend):
        contig = str(xstart // 1000000000)
        start = xstart % 1000000000
        end = xend % 1000000000
        return self.getCoverage(contig, start, end) 


class CoverageCollection(object):

    def __init__(self):
        self.collection = IntervalTree()      

    def setTabixPath(self, min_length, max_length, contig, path):
        coverages = self.collection.search(min_length, max_length) 
        if coverages:
            for coverage in coverages:
                coverage.data.setTabixPath(contig, path)
                break
        else:
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
            return coverage.data.getCoverage(contig, start, end)

        return []

    def getCoverageX(self, xstart, xend):
        for coverage in self.collection.search(xend - xstart):
            return coverage.data.getCoverageX(xstart, xend)

        return []     

    def printAll(self):
        for coverage in self.collection.items():
            print coverage
