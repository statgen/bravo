#!/bin/bash

# Build tools and load data into MongoDB. Insert proper container name
docker exec bravo_web_1 /bin/bash -c /deploy/INSTALL.sh -v hg38 -t 8

# Create variants collection in MongoDB.
docker exec bravo_web_1 python manage.py variants -t 8 -v /data/import_vcf/chr22.TOPMed_freeze5_62784.vcf.gz
