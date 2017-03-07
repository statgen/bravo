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


### Make some config for your new dataset
In `flask_config.py`, there's one section (actually a python class) for each dataset.  Make (or repurpose) a section for your dataset.

Then make sure that `exac.py` uses the name of your dataset on the line `app.config.from_object('flask_config.<name_of_dataset>')`.


## Set up OAuth and an email whitelist
In your section of `flask_config.py`, the variable `EMAIL_WHITELIST` should be a list of allowed email addresses.  Currently that list is made in a separate file like `whitelist_topmed.py` and imported into `flask_config.py`, but you could just use a list instead.  If the list is empty or false (ie, `EMAIL_WHITELIST = False`), any email will be allowed.

You need to set up a OAuth with Google.  Go [here](https://console.developers.google.com/apis/credentials) and create a project.  In the list "Authorized redirect URIs" add your OAuth callback URL, which should look like `https://example.com/callback/google` or `https://example.com:5000/callback/google`.  Then copy the client ID and secret from the top of that page into `flask_config.py` for the variables `GOOGLE_LOGIN_CLIENT_ID` and `GOOGLE_LOGIN_CLIENT_SECRET`.


## Set up Google Analytics (optional)
Go [here](https://analytics.google.com/analytics/web) and do whatever you have to to get your own `UA-xxxxxx-xx` tracking id.  Put that in `flask_config.py`.  Or just leave the default `UA-01234567-89`, and you won't receive any of the tracking data.


### Import data into Mongo

You'll need the following files:

- `ALL.vcf.gz` and `ALL.vcf.gz.tbi`
    - this is the VCF from the section "Prepare a VCF"
    - stored in the variable `SITES_VCFS`.  If you want to use a different name or multiple files, just make sure that pattern for `SITES_VCFS` in `flask_config.py` matches all of your files.
    - used by `load_variants_file()`

- `canonical_transcripts.txt.gz`
    - stored in the variable `CANONICAL_TRANSCRIPT_FILE`
    - used by `load_gene_models()`
    - I'm not sure where to get this.  I see that exac provides one [here](https://personal.broadinstitute.org/konradk/exac_browser/), but that's old.  I can't understand [this page](https://www.gencodegenes.org/gencode_tags.html).  Gencode talks about canonical transcripts [here](http://www.ensembl.org/Help/Glossary?id=346;redirect=no) but I can't find any other references to them.  APPRIS lists canonical transcripts [here](http://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/GRCh37/appris_data.principal.txt) but some genes are missing and some have two tied principle transcripts.  Maybe I just need to follow Gencode's definition and make these by hand.

- `dbsnp149.txt.bgz` and `dbsnp149.txt.bgz.tbi`
    - stored in `DBSNP_FILE`
    - used by `load_dbsnp_file()`
    - as dbSNP updates, they drop older versions, so you might have to update from 149 to something higher. Everything should probably work fine with newer versions of dbSNP
    - you should be able to make this by running (possibly after updating version numbers):

            wget -O dbsnp149.bcp.gz ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b149_GRCh37p13/database/organism_data/b149_SNPChrPosOnRef_105.bcp.gz
            zcat dbsnp149.bcp.gz | awk '$3 != ""' | perl -pe 's/ +/\t/g' | sort -k2,2 -k3,3n | bgzip -c > dbsnp149.txt.bgz
            tabix -s 2 -b 3 -e 3 dbsnp149.txt.bgz

    - you will have to update the pattern for this file in your section of `flask_config.py`.

- `gencode.gtf.gz`
    - stored in `GENCODE_GTF`
    - used by `load_gene_models()`
    - this can be downloaded from [here](https://personal.broadinstitute.org/konradk/exac_browser/), but that's out-of-date.  Gencode v25 for GRCh37 is [here](ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gtf.gz) and should work fine but hasn't been tested.

- `omim_info.txt.gz`
    - stored in `OMIM_FILE`
    - used by `load_gene_models()`
    - This can be downloaded from [ExAC](https://personal.broadinstitute.org/konradk/exac_browser/), but that's out-of-date.  This isn't [mim2gene](https://omim.org/static/omim/data/mim2gene.txt), maybe it's one of the educational-use-only files on that download page.

- `dbNSFP2.6_gene.gz`
    - stored in `DBNSFP_FILE`
    - used by `load_gene_models()`
    - this can be downloaded from [this site](https://sites.google.com/site/jpopgen/dbNSFP).  Our version is out-of-date and should be replaced by <http://genenames.org>.

Put these in a directory, and store that directory in `_FILES_DIRECTORY` in your section of `flask_config.py`.  If you use exactly these names, everything should work.  If you change names, modify the file-matching patterns in your section of `flask_config.py`.

Then run:

    ./manage.py load_variants_file
    ./manage.py load_dbsnp_file
    ./manage.py load_gene_models


### Precalculate some things in the database

    ./manage.py precalculate_metrics
    ./manage.py precalculate_whether_variant_is_ever_missense_or_lof


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
    ./exac.py

And visit on your browser:

    http://localhost:5000
    http://localhost:5000/gene/ENSG00000237683
    http://localhost:5000/variant/20-76735-A-T

For testing, you can open up an interactive shell with:

    ./manage.py shell
