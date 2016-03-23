Installation
============

### Dependencies

Install MongoDB.  Configure it and make it run all the time.

Install python packages in a virtualenv.

    virtualenv venv
    source venv/bin/activate
    pip2 install -r requirements.txt

Some packages will require Python headers (python-dev on some systems).

You probably want to run Flask behind something else, like Apache2.


### Prepare a VCF

Use the scripts in `data/new_vcf_info/` to extract a subset of samples from a vcf file and calculate a few summaries for each variant.

Use `data/add_cadd_scores.py` to add CADD scores to a vcf.

Use `data/remove_ac0.py` to remove variants that never vary.

Use `data/import_info.py` with some options (I don't know which) to copy some INFO fields from one vcf (eg, your full vcf) into another (eg, your sites vcf)


### Import data into Mongo

You'll need files like the following:

    ALL.polymorphic.vcf.gz (this is the VCF from the section "Prepare a VCF")
    ALL.polymorphic.vcf.gz.tbi
    canonical_transcripts.txt.gz
    dbsnp144.txt.bgz
    dbsnp144.txt.bgz.tbi
    gencode.gtf.gz
    omim_info.txt.gz
    dbNSFP2.6_gene.gz (used in parsing.py? TODO)

Put these in a directory, and store that directory in `_FILES_DIRECTORY` in `flask_config.py`.  As long as the names match our patterns, those files should correctly be stored in:

- `SITES_VCFS`: used by `load_variants_file()`
- `DBSNP_FILE`: used by `load_dbsnp_file()`
- `GENCODE_GTF`: used by `load_gene_models()`
- `CANONICAL_TRANSCRIPT_FILE`: used by used by `load_gene_models()`
- `OMIM_FILE`: used by `load_gene_models()`

Then run:

    ./manage.py load_variants_file
    ./manage.py load_dbsnp_file
    ./manage.py load_gene_models


### Precalculate some things in the database

    python manage.py precalculate_metrics
    python manage.py precalculate_whether_variant_is_ever_missense_or_lof


### Generate coverage files and configure the browser to use them

1. Use the code in `data/base_coverage/glf2depth/` to create a full coverage file (ie, with coverage for every available base).
Make one `*.full.json.gz` for each chromosome in some directory.

2. Use the scripts in `data/base_coverage/` to bin the coverage.
Make a couple directories with different levels of binning (and again, one `.json.gz` file for each chromosome).

3. Tabix them all.

4. Reference all of the coverage files in `BASE_COVERAGE` in `flask_config.py`.


### Create the user table:

While the other operations here are all idempotent, this one will wipe your user data, so only run it when you don't yet have user data.

    ./manage.py create_users


### Start the server

You can run the development server with:

    source venv/bin/activate
    python2 exac.py

And visit on your browser:

    http://localhost:5000
    http://localhost:5000/gene/ENSG00000237683
    http://localhost:5000/variant/20-76735-A-T

For testing, you can open up an interactive shell with:

    python manage.py shell
