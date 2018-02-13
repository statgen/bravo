import os
import imp
import inspect
import warnings
import pymongo
import re
import time

def load_config(name):
    f, path, desc = imp.find_module(name.strip().split('.')[0])
    m = imp.load_module(name, f, path, desc)
    config = dict()
    for member, value in inspect.getmembers(getattr(m, name.strip().split('.')[1])):
        config[member] = value
    return config

def delete_cache(db, collection, cursor, delay):
    ids_to_remove = []
    files_to_remove = []
    for l in cursor:
        ids_to_remove.append(l['_id'])
        files_to_remove.append((l['bam'], l['bai']))
    if len(ids_to_remove) > 0:
        requests = [ pymongo.DeleteMany({ '_id': { '$in': ids_to_remove } })]
        result = db[collection].bulk_write(requests)
    if len(files_to_remove) > 0: 
        time.sleep(delay) # sleep a little to allow finish any active downloads of selected cached files          
        for bam, bai in files_to_remove:
            try:
                os.remove(os.path.join(cache_dir, bam))
                os.remove(os.path.join(cache_dir, bai))
            except OSError:
                pass

if __name__ == '__main__':
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        config = load_config('flask_config.BravoFreeze5GRCh38Config')
    mongo_host = config['MONGO']['host']
    mongo_port = config['MONGO']['port']
    mongo_db_name = config['MONGO']['name']
    cache_collection = config['IGV_CACHE_COLLECTION']
    cache_dir = config['IGV_CACHE_DIRECTORY']
    cache_limit = config['IGV_CACHE_LIMIT']

    mongo = pymongo.MongoClient(host = mongo_host, port = mongo_port, connect = False)
    db = mongo[mongo_db_name]

    cached_files_filter = re.compile('(het|hom)-', re.IGNORECASE)
    deleted_cached_files_filter = re.compile('delete-', re.IGNORECASE)

    n = db[cache_collection].count({ 'name': cached_files_filter })

    if n >= cache_limit:
        n_to_remove = int(0.5 * cache_limit) + (n - cache_limit)
        to_remove = db[cache_collection].find({ 'name': cached_files_filter }, sort = [('accessed', pymongo.ASCENDING)], limit = n_to_remove)
        delete_cache(db, cache_collection, to_remove, 15)

    to_remove = db[cache_collection].find({ 'name': deleted_cached_files_filter })
    delete_cache(db, cache_collection, to_remove, 0)
