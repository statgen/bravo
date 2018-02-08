import os
import pysam
import pymongo
from utils import Xpos

class SequencesClient(object):
    '''Manages CRAMS for all chromosomes. Assumes one CRAM per chromosome.'''
    def __init__(self, crams_dir, reference_path, window_bp):
        self._crams_dir = crams_dir
        self._reference_path = reference_path
        self._window_bp = window_bp
        self._crams = dict()
        for cram_file in os.listdir(self._crams_dir):
            if cram_file.endswith(".cram"):
                cram_path = os.path.join(self._crams_dir, cram_file)
                with pysam.AlignmentFile(cram_path, 'rc', reference_filename = self._reference_path) as icram:
                    chrom = None
                    for read in icram:
                        chrom = read.reference_name
                        break
                    self._crams[chrom] = { 'header': { 'HD': icram.header['HD'], 'SQ': icram.header['SQ'] }, 'path': cram_path }
    def get_samples(self, db, variant_id):
        chrom, pos, ref, alt = variant_id.split('-')
        pos = int(pos)
        xpos = Xpos.from_chrom_pos(chrom, pos)
        variant = db.variants.find_one({'xpos': xpos, 'ref': ref, 'alt': alt}, projection = {'_id': False})
        if variant is None:
            return None
        cram = self._crams.get(chrom, None)
        if cram is None:
            chrom = chrom[3:] if chrom.startswith('chr') else 'chr{}'.format(chrom)
            cram = self._crams.get(chrom, None)
            if cram is None:
                return None
        start = pos - self._window_bp if pos > self._window_bp else 0
        stop = pos + self._window_bp
        qname = '{}:{}:{}:'.format(pos, ref, alt)
        samples = set()
        with pysam.AlignmentFile(cram['path'], 'rc', reference_filename = self._reference_path) as icram:
            for read in icram.fetch(chrom, start, stop):
                if read.query_name.startswith(qname):
                    sample = read.query_name.split(':')[3]
                    samples.add('{}-{}'.format('hom' if sample.startswith('0') else 'het', int(sample)))
        return { 'names': sorted(samples) }
    def get_bai(self, db, variant_id, sample_id):
        chrom, pos, ref, alt = variant_id.split('-')
        pos = int(pos)
        xpos = Xpos.from_chrom_pos(chrom, pos)
        variant = db.variants.find_one({'xpos': xpos, 'ref': ref, 'alt': alt}, projection = {'_id': False})
        if variant is None:
            return None
        cram = self._crams.get(chrom, None)
        if cram is None:
            chrom = chrom[3:] if chrom.startswith('chr') else 'chr{}'.format(chrom)
            cram = self._crams.get(chrom, None)
            if cram is None:
                return None
        start = pos - self._window_bp if pos > self._window_bp else 0
        stop = pos + self._window_bp
        sample_type, sample_no = sample_id.split('-')
        qname = '{}:{}:{}:{}{}:'.format(pos, ref, alt, 0 if sample_type == 'hom' else '', sample_no)
        filename = '{}.{}.bam'.format(variant_id, sample_id)
        print qname, filename
        with pysam.AlignmentFile(cram['path'], 'rc', reference_filename = self._reference_path) as icram, pysam.AlignmentFile(filename, 'wb', header = cram['header']) as obam:
            for read in icram.fetch(chrom, start, stop):
                if read.query_name.startswith(qname):
                    for tag, value in read.get_tags():
                        read.set_tag(tag, None)
                    obam.write(read)
        return filename
    def get_bam(self, db, variant_id, sample_id, start, end = None):
        filename = '{}.{}.bam'.format(variant_id, sample_id)
        file_size = os.path.getsize(filename)
        start = int(start)
        length = int(end or file_size) - start
        result = { 'data': None, 'start': start, 'end': start + length - 1, 'size': file_size }
        with open(filename, 'rb') as f:
            f.seek(start)
            result['data'] = f.read(length)
        return result
        
