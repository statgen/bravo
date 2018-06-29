#!/bin/bash

mongodump -d <database> -c <collection> --archive  --gzip | gsutil cp - gs://<bucket>/test.bson.gz
