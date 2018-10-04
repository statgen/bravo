#!/bin/bash

# Dump backup from mongo container and send to google cloud storage bucket
docker exec -i <container-name> sh -c 'exec mongodump -d <database> --archive --gzip' | gsutil cp - gs://<bucket>/test.archive.gz

# Restore from backup file after downloading from bucket
docker exec -i <container-name> sh -c 'exec mongorestore -d <database> --archive --gzip' < /data/test.archive.gz
