#!/usr/bin/env python2

"""
Read files and create `db.substudies` for TOPMed.

Abilities:
1. Map variant -> NWD_IDs (/var/browser_genotypes/topmed_freeze2/topmed_freeze2_10597.chr{chrom}.overlap_removed.svm_pass.genotypes.vcf.gz)
2. Map NWD_ID -> study (/var/browser_genotypes/topmed_freeze2/lookup.table.tab) (TODO is this useful?)
3. Map NWD_ID -> sub-study and sub-study internal ID (/net/topmed/incoming/study.reference/phs)
4. Map sub-study -> users who are contacted and may accept/decline (PI_email_addresses.txt)

Database:
db.substudies = [{substudy_name, {nwd_id:internal_id, ...}}, ...]

# TODO
- [ ] Is there still any value in /var/browser_genotypes/topmed_freeze2/lookup.table.tab?  Check the overlap of those studies with the other studies. Then ask Tom.

To watch for new files, use:
while diff -q <(ls | grep -o '^[0-9]*' | sort -u) ~/tmp/dates.txt; do sleep 600; done; echo | mail -s 'new files in /net/topmed/incoming/study.reference/phs' pjvh@umich.edu ; ls | grep -o '^[0-9]*' | sort -u > ~/tmp/dates.txt
"""

import glob
import gzip
import xml.etree.ElementTree as ET

import pymongo

import os.path
import sys
sys.path.insert(0, os.path.abspath(os.path.join(__file__, '../../..')))
import flask_config

def get_db():
    MONGO_HOST = flask_config.BravoConfig.MONGO['host']
    MONGO_PORT = flask_config.BravoConfig.MONGO['port']
    MONGO_COLLECTION_NAME = flask_config.BravoConfig.MONGO['name']

    mongo_client = pymongo.MongoClient(host=MONGO_HOST, port=MONGO_PORT)
    return mongo_client[MONGO_COLLECTION_NAME]
db = get_db()
db.substudies.drop()

def get_filenames():
    most_recent_date = max(os.path.basename(s)[:8] for s in glob.glob('/net/topmed/incoming/study.reference/phs/20*'))
    assert most_recent_date.isdigit(), most_recent_date
    recent_files = glob.glob('/net/topmed/incoming/study.reference/phs/{}*'.format(most_recent_date))
    assert 10 < len(recent_files) < 100
    assert all(filename.endswith('.xml.gz') for filename in recent_files)
    return recent_files
for fname in get_filenames():
    # TODO: parse PI_name from fname
    # print(fname)
    with gzip.open(fname) as f:
        root = ET.fromstring(f.read())

        assert len(root.findall('.//Study')) == 1
        substudy_elem = root.find('.//Study')
        substudy = {
            'substudy_name': substudy_elem.get('study_handle'),
            'substudy_longname': substudy_elem.get('study_name'),
            'samples': [],
        }
        for sample_elem in root.findall('.//Sample'):
            if sample_elem.get('submitted_sample_id').startswith('SM-'): continue # I don't know what these are.
            # assert sample_elem.get('repository').lstrip('NHLBI_').lstrip('TOPMed_WGS_') == substudy['substudy_name'].lstrip('NHLBI_').lstrip('TOPMed_WGS_'), \
            #     (sample_elem.get('repository'), substudy['substudy_name'])
            assert sample_elem.get('submitted_sample_id').startswith('NWD'), sample_elem.get('submitted_sample_id')
            assert sample_elem.get('sex') in ['male', 'female', None], sample_elem.get('sex')
            sample = {
                'NWD_ID': sample_elem.get('submitted_sample_id'),
                'submitted_subject_id': sample_elem.get('submitted_subject_id'),
                'consent_short_name': sample_elem.get('consent_short_name'),
            }
            sex = sample_elem.get('sex')
            if sex is not None:
                sample['sex'] = sex

            substudy['samples'].append(sample)

    db.substudies.insert(substudy)
