Installation
============

1. [System Dependencies](#system-dependencies)
2. [Data Preparation](#data-preparation)
   1. [Prepare VCF](#prepare-vcf)
   2. [Prepare percentiles](#prepare-percentiles)
   3. [Prepare coverage](#prepare-coverage)
   4. [Prepare CRAM](#prepare-cram)
3. [Data Import](#data-import)
4. [Authentication](#authentication)
5. [Email Whitelist](#email-whitelist)
6. [Google Analytics](#google-analytics)
7. [Start the server](#start-the-server)

## System Dependencies

Install MongoDB.  Configure it and make it run all the time.

Install python packages in a virtualenv.

    virtualenv venv
    source venv/bin/activate
    pip2 install -r requirements.txt

Some packages will require Python headers (python-dev on some systems).

You probably want to run Flask behind something else, like Apache2.

## Data Preparation

In `data/` directory you will find tools/scripts to prepare your data for importing into Mongo database and using in Bravo browser.
Some of these tools are implemented in C++ and need to be compiled as follows:

    mkdir data/DataPrep/build
    cd data/DataPrep/build
    cmake ..
    make install

After successfull compilation, the executables will be installed in `data/DataPrep/build/bin`.
Table below lists all tools that are needed for data preparation.

| Tool | Localtion | Description |
|:-----|:----------|:------------|
| CountAlleles | `data/DataPrep/build/bin` | Computes NS, AN, AC, AF, Hom, Het values for each variant in a subset of samples (monomorphic variants are dropped) |
| GTHistogram | `data/DataPrep/build/bin` | Computes histograms for DP and GQ fields for each variant in a subset of samples |
| INFOPercentiles | `data/DataPrep/build/bin` | This program computes percentiles for the QUAL field or any arbitrary numeric INFO field across all variants |
| prepare_sequences.py | `data` | Generates CRAM file with sequences from heterozygous/homozygous samples |

Specify `--help` to see all possible options for each of the above tools (e.g. `Count Alleles --help` or `python prepare_sequences.py --help`).

In addition to the in-house tools listed above, you will need Variant Effect Predictor (VEP) and bcftools.

### Prepare VCF

If applicable, prepare a list of sample ID's that you want to include into Bravo.
Then, on your raw VCF file(s) run the following (in parallel):

1. Variant Effect Predictor (VEP). Genotypes are not required.
2. CountAlleles. Requires genotypes (GT format field).
3. GTHistograms. Requires genotypes (GT format field), depth (DP format field), and quality (GQ format field).

Each tool will produce VCF file(s) with computed INFO fields. At the final step, merge this info fields into a single VCF using `bcftools annotate` command.

### Prepare percentiles

To compute percentiles for each INFO field in each variant, use `INFOPercentiles` tool.

### Prepare coverage

1. Use the code in `data/base_coverage/glf2depth/` to create a full coverage file (ie, with coverage for every available base).
Make one `*.full.json.gz` for each chromosome in some directory.

2. Use the scripts in `data/base_coverage/` to bin the coverage.
Make a couple directories with different levels of binning (and again, one `.json.gz` file for each chromosome).

3. Tabix them all.

4. Reference all of the coverage files in `BASE_COVERAGE` in `flask_config.py`.

### Prepare CRAM

Bravo uses IGV.js to visualize raw sequenced from up to 10 random alternate allele carriers.

## Data Import

### Make some config for your new dataset
In `flask_config.py`, there's one section (actually a python class) for each dataset.  Make (or repurpose) a section for your dataset.

Then make sure that `exac.py` uses the name of your dataset on the line `app.config.from_object('flask_config.<name_of_dataset>')`.

### Import data into Mongo

You'll need the following files:

- `ALL.vcf.gz` and `ALL.vcf.gz.tbi`
    - this is the VCF from the section "Prepare a VCF"
    - stored in the variable `SITES_VCFS`.  If you want to use a different name or multiple files, just make sure that pattern for `SITES_VCFS` in `flask_config.py` matches all of your files.
    - used by `load_variants_file()`

- `canonical_transcripts.txt.gz`
    - stored in the variable `CANONICAL_TRANSCRIPT_FILE`
    - used by `load_gene_models()`
    - this can be obtained with `perl download_canonical_transcripts.pl | gzip -c > canonical_transcripts.txt.gz`. Note that to run `download_canonical_transcripts.pl` you will need to install Ensembl API perl modules.

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
    - this can be downloaded from [here](https://www.gencodegenes.org/releases/current.html)

- `omim_info.txt.gz`
    - stored in `OMIM_FILE`
    - used by `load_gene_models()`
    - this can be downloaded from Ensembl BioMart:
        - "Ensembl Genes" database 
        - "Human genes" dataset 
        - "Gene stable ID", "Transcript stable ID", "MIM gene accession", "MIM gene description" attributes

- `dbNSFP2.6_gene.gz`
    - stored in `DBNSFP_FILE`
    - used by `load_gene_models()`
    - this can be downloaded from [this site](https://sites.google.com/site/jpopgen/dbNSFP).  Our version is out-of-date and should be replaced by <http://genenames.org>.

Put these in a directory, and store that directory in `_FILES_DIRECTORY` in your section of `flask_config.py`.  If you use exactly these names, everything should work.  If you change names, modify the file-matching patterns in your section of `flask_config.py`.

Then run:

    ./manage.py load_variants_file
    ./manage.py load_dbsnp_file
    ./manage.py load_gene_models


### Create the user table:

While the other operations here are all idempotent, this one will wipe your user data, so only run it when you don't yet have user data.

    ./manage.py create_users
    
## Authentication

Bravo supports user authentication using Google's OAuth 2.0 protocol, which is optional and is disabled by default. This section describes how to enable it.

First, make sure that your Bravo instance is served using HTTPS protocol.

Second, you need to set up a OAuth with Google. Go [here](https://console.developers.google.com/apis/credentials) and create a project. Your project will get a `Client ID` and a `Client secret`. In the list "Authorized redirect URIs" add your OAuth callback URL, which should look like `https://[bravo base URL]/callback/google` (e.g. `https://mybravo.myinstitution.org/callback/google`).

**Attention!** Don't expose to anyone your `Client ID` and `Client secret`, and make sure you are using HTTPS for your callback URL.

Third, follow these steps to enable authentication in Bravo:
1. Set the `GOOGLE_AUTH` variable in Bravo configuration file to `True`.
2. Assign the `GOOGLE_LOGIN_CLIENT_ID` variable in Bravo configuration file to your `Client ID` from Google.
3. Assign the `GOOGLE_LOGIN_CLIENT_SECRET` variable in Bravo configuration file to your `Client secret` from Google.


## Email Whitelist

Bravo allows whitelist specific users based on their email address. To enable whitelisting, follow these steps:
1. Set up user authentication as described in [Authentication](#authentication) section.
2. Set the `EMAIL_WHITELIST` variable in Bravo configuration file to `True`.
3. Import list of emails from a text file (one email per line) to the Mongo database:
    
       ./manage.py load_whitelist -w emails.txt
       

## Google Analytics
This step is optional. Go [here](https://analytics.google.com/analytics/web) and do whatever you have to to get your own `UA-xxxxxx-xx` tracking id.  Put that in `flask_config.py`.  Or just leave the default `UA-01234567-89`, and you won't receive any of the tracking data.

## Start the server

You can run the development server with:

    source venv/bin/activate
    ./exac.py

And visit on your browser:

    http://localhost:5000
    http://localhost:5000/gene/ENSG00000237683
    http://localhost:5000/variant/20-76735-A-T

For testing, you can open up an interactive shell with:

    ./manage.py shell
