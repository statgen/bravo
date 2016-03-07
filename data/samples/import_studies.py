#!/usr/bin/env python2

"""
Read files and create `db.studies` for TOPMed.
"""

import pymongo

import os.path
import sys
sys.path.insert(0, os.path.abspath(os.path.join(__file__, '../..')))
import flask_config

MONGO_HOST = flask_config.BravoConfig.MONGO['host']
MONGO_PORT = flask_config.BravoConfig.MONGO['port']
MONGO_COLLECTION_NAME = flask_config.BravoConfig.MONGO['name']

mongo_client = pymongo.MongoClient(host=MONGO_HOST, port=MONGO_PORT)
db = mongo_client[MONGO_COLLECTION_NAME]
db.studies.drop()

with open('PI_email_addresses.txt') as f:
    for line in f:
        pi_name, study_name, member_emails = line.split()
        if member_emails == "HAPMAP":
            db.studies.insert({
                'study_name': study_name,
                'HAPMAP': True,
            })
        else:
            member_emails = member_emails.split(',')
            assert all('@' in email for email in member_emails)
            db.studies.insert({
                'study_name': study_name,
                'member_emails': member_emails,
            })
