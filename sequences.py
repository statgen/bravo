import os
import random
import string

import pymongo
import pysam
from utils import Xpos


class SequencesClient(object):
    '''Manages CRAMS for all chromosomes. Assumes one CRAM per chromosome.'''
    def __init__(self, crams_dir, reference_path, cache_dir, cache_collection, window_bp):
        self._crams_dir = crams_dir
        self._reference_path = reference_path
        self._window_bp = window_bp
        self._crams = dict()
        if not os.path.exists(cache_dir):
            raise Exception('Provided cache path does not exist.')
        if not os.path.isdir(cache_dir):
            raise Exception('Provided cache path must be a directory.')
        try: # check if has write permissions in cache directory
            filename = os.path.join(cache_dir, SequencesClient.get_random_filename(10))
            f = open(filename, 'w')
        except:
            raise Exception('Error while writing to the provided cache directory.')
        else:
            f.close()
            os.remove(filename)
        self._cache_dir = cache_dir
        if len(cache_collection.strip()) == 0:
            raise Exception('Cache collection name cannot be empty.')
        self._cache_collection = cache_collection
        for cram_file in os.listdir(self._crams_dir):
            if cram_file.endswith('.cram'):
                cram_path = os.path.join(self._crams_dir, cram_file)
                with pysam.AlignmentFile(cram_path, 'rc', reference_filename = self._reference_path) as icram:
                    chrom = None
                    for read in icram:
                        chrom = read.reference_name
                        break
                    self._crams[chrom] = { 'header': { 'HD': icram.header['HD'], 'SQ': icram.header['SQ'] }, 'path': cram_path }

    @staticmethod
    def get_random_filename(length):
        return ''.join(random.choice(string.ascii_letters + string.digits) for x in xrange(length))

    @staticmethod
    def create_cache_collection_and_index(db, collection_name):
        if collection_name in db.collection_names():
            db[collection_name].drop()
        db.create_collection(collection_name)
        db[collection_name].create_index('name', unique = True)

    def create_bam(self, db, variant_id, sample_id):
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
        # check if exists in cache
        cache_name = '{}.{}'.format(variant_id, sample_id)
        cache_entry = db[self._cache_collection].find_one({ 'name': cache_name }, projection = { '_id': False })
        if cache_entry is not None:
            if os.path.exists(cache_entry['bam']) and os.path.exists(cache_entry['bai']):  # check if files still exist in cache directory
                db[self._cache_collection].update_one({ 'name': cache_name }, { '$inc': { 'accessed' : 1 } })
                return { 'bam': cache_entry['bam'], 'bai': cache_entry['bai'] }
            else:
                delete_cache_name = 'delete-{}'.format(SequencesClient.get_random_filename(10))
                db[self._cache_collection].update_one({ 'name': cache_name }, { '$set': { 'name': delete_cache_name } })
        start = pos - self._window_bp if pos > self._window_bp else 0
        stop = pos + self._window_bp
        sample_type, sample_no = sample_id.split('-')
        bam_path = os.path.join(self._cache_dir, '{}.{}.{}.bam'.format(variant_id, sample_id, SequencesClient.get_random_filename(5)))
        bai_path= '{}.bai'.format(bam_path)
        qname = '{}:{}:{}:{}{}:'.format(pos, ref, alt, 0 if sample_type == 'hom' else '', sample_no)
        with pysam.AlignmentFile(cram['path'], 'rc', reference_filename = self._reference_path) as icram, pysam.AlignmentFile(bam_path, 'wb', header = cram['header']) as obam:
            for read in icram.fetch(chrom, start, stop):
                if read.query_name.startswith(qname):
                    for tag, value in read.get_tags():
                        read.set_tag(tag, None)
                    obam.write(read)
        pysam.index(bam_path)
        # save to cache if does not exist yet
        try:
            result = db[self._cache_collection].update_one({ 'name': cache_name}, { '$inc': { 'accessed': 1 }, '$setOnInsert': { 'bam': bam_path, 'bai': bai_path } }, upsert = True)
        except pymongo.errors.DuplicateKeyError:
            delete_cache_name = 'delete-{}'.format(SequencesClient.get_random_filename(10))
            db[self._cache_collection].insert_one({'name': delete_cache_name, 'bam': bam_path, 'bai': bai_path })
        return { 'bam': bam_path, 'bai': bai_path }

    def get_samples(self, db, variant_id):
        chrom, pos, ref, alt = variant_id.split('-')
        pos = int(pos)
        xpos = Xpos.from_chrom_pos(chrom, pos)
        variant = db.variants.find_one({'xpos': xpos, 'ref': ref, 'alt': alt}, projection = {'_id': False})
        if variant is None:
            return None
        # check if exists in cache
        cache_entry =  db[self._cache_collection].find_one({'name': variant_id}, projection = {'_id': False})
        if cache_entry is not None:
            return { 'names': cache_entry['sample'] }
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
        samples = sorted(samples, key = lambda x: int(x.split('-')[1]))
        # save to cache if does not exist yet
        try:
            result = db[self._cache_collection].update_one({'name': variant_id},  { '$setOnInsert': { 'sample': samples } }, upsert = True)
        except pymongo.errors.DuplicateKeyError:
            pass
        return { 'names': samples }

    def get_bai(self, db, variant_id, sample_id):
        bam = self.create_bam(db, variant_id, sample_id)
        return bam['bai'] if bam else None

    def get_bam(self, db, variant_id, sample_id, start, end = None):
        bam = self.create_bam(db, variant_id, sample_id)
        if not bam:
            return None
        file_size = os.path.getsize(bam['bam'])
        start = int(start)
        length = int(end or file_size) - start
        result = { 'data': None, 'start': start, 'end': start + length - 1, 'size': file_size }
        with open(bam['bam'], 'rb') as f:
            f.seek(start)
            result['data'] = f.read(length)
        return result
