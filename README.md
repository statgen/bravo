# BRAVO
BRowse All Variants Online

## Installation

1. [System Setup](#system-setup)
1. [Configuration](#configuration)
    - [Application Settings](#application-Settings)
    - [Apache Config](#apache-config)
    - [Access Control](#access-control)
        - [Authentication](#authentication)
        - [Email Whitelist](#email-whitelist)
        - [Terms of Use](#terms-of-use)
    - [Google Analytics](#google-analytics)
1. [Launch the Application](#launch-the-application)
1. [Data Preparation](#data-preparation)
    1. [Prepare VCF](#prepare-vcf)
    1. [Prepare percentiles](#prepare-percentiles)
    1. [Prepare coverage](#prepare-coverage)
    1. [Prepare CRAM](#prepare-cram)
1. [Load Data]()
1. [Data Backup and Restore](#data-backup-and-restore)

## System Setup

BRAVO is packaged using [Docker](https://www.docker.com/). This is the recommended way by which to deploy the application.
System dependencies can be installed using the `setup.sh` script. The script installs and configures the following applications on your machine:

- [Docker CE](https://www.docker.com/products/docker-engine)
- [Docker Compose](https://docs.docker.com/compose/)
- [Google Cloud SDK](https://cloud.google.com/sdk/)(Used for application data backups.)

The script will also create the default directory structure on the local host at `/data`.

```bash
|-- data
|   |-- cache
|   |   |-- igv_cache
|   |-- coverage
|   |-- cram
|   |-- genomes
|   |-- import_vcf
```

## Configuration

### Application Settings 

The BRAVO configuration file for project dependent settings can be found at `config/default.py`.
You will need to edit some of the settings to match your environment.

### Apache Config

It is recommended that you run your BRAVO application behind a reverse proxy such as Apache. An example Apache configuration can be found in `apache-example.conf`.
In order to use Apache as a reverse proxy you will need to enable the `mod_proxy` and `mod_proxy_http` modules using the command `sudo a2enmod proxy proxy_http`.

You will also need to change the BRAVO configuration file setting `PROXY` from `False` to `True` and reload your application as described in the [Launch the Application](#launch-the-application) section.

### Access Control

#### Authentication

BRAVO supports user authentication using Google's OAuth 2.0 protocol, which is optional and is disabled by default. This section describes how to enable it.

First, make sure that your BRAVO instance is served using HTTPS protocol.

Second, you need to set up a OAuth with Google. Go [here](https://console.developers.google.com/apis/credentials) and create a project. Your project will get a `Client ID` and a `Client secret`. In the list "Authorized redirect URIs" add your OAuth callback URL, which should look like `https://[bravo base URL]/callback/google` (e.g. `https://mybravo.myinstitution.org/callback/google`).

**Attention!** Don't expose to anyone your `Client ID` and `Client secret`, and make sure you are using HTTPS for your callback URL.

Third, follow these steps to enable authentication in BRAVO: 
1. Set the `GOOGLE_AUTH` variable in the  BRAVO configuration file to `True`.
2. Assign the `GOOGLE_LOGIN_CLIENT_ID` variable in the BRAVO configuration file to your `Client ID` from Google.
3. Assign the `GOOGLE_LOGIN_CLIENT_SECRET` variable in the BRAVO configuration file to your `Client secret` from Google.

#### Email Whitelist

BRAVO allows whitelist specific users based on their email address. To enable whitelisting, follow these steps:
1. Set up user authentication as described in [Authentication](#authentication) section.
2. Set the `EMAIL_WHITELIST` variable in BRAVO configuration file to `True`.
3. Import list of emails from a text file (one email per line) to the Mongo database:
    
       ./manage.py whitelist -w emails.txt
 
#### Terms of Use

If your BRAVO users must to agree to any terms/conditions before browsing your data, you need to enable **Terms of Use** page as follows:
1. Set up user authentication as described in [Authentication](#authentication) section.
2. Set the `TERMS` variable in BRAVO configuration file to `True`.
3. Write your terms/conditions to the `templates/terms.html` file.


### Google Analytics

This step is optional. Go [here](https://analytics.google.com/analytics/web) and do whatever you have to to get your own `UA-xxxxxx-xx` tracking id.  Put that in `default.py`.  Or just leave the default `UA-01234567-89`, and you won't receive any of the tracking data.

## Launch the Application

In order to start the service, run `docker-compose up -d` from the projects home directory.

In order to stop the application run `docker-compose down` from the projects home directory.

If you are making changed to configuration files or code and need to reload the changes, run `docker-compose down && docker-compose up --build -d` in order to stop rebuild the containers with your updates.

## Data Preparation

In the `data/` directory you will find tools/scripts to prepare your data for importing into Mongo database and using in BRAVO browser.

1. Compile data preparation tools:
   ```
   cd data/DataPrep/
   cget install.
   ```
   After successful compilation, the executables will be installed in `data/DataPrep/cget/bin`.
   
2. Preprare VCF with the INFO fields NS, AN, AC, AF, Hom, Het, DP, AVGDP, AVGDP_R, AVGGQ, AVGGQ_R, DP_HIST, DP_HIST_R, GQ_HIST, and GQ_HIST_R:
   ```
   ./ComputeAlleleCountsAndHistograms -i [input bcf/vcf] -s [samples file] -r [CHR:START-END] -o [output.vcf.gz]
   ```
   Input BCF/VCF must have DP and GQ FORMAT fields. Input BCF/VCF can be accessed both from local and from Google bucket storages. The input samples file (one sample ID per line) and chromosomal region CHR:START-END are optional.

3. Run [Variant Effect Predictor (VEP)](https://www.ensembl.org/vep) on the VCF created in step (2):
   ```
   ./vep -i [input vcf.gz] --plugin LoF --assembly [GRCh37/GRCh38] --cache --offline --vcf --sift b --polyphen b --ccds --uniprot --hgvs --symbol --numbers --domains --regulatory --canonical --protein --biotype --af --af_1kg --pubmed --shift_hgvs 0 --allele_number --format vcf --force --buffer_size 100000 --compress_output gzip --no_stats -o [output vcf.gz]
   ```
   Note: specify [LoF plugin](https://github.com/konradjk/loftee) configuration options as you need.

Table below lists all tools that are needed for data preparation.

| Tool | Location | Description |
|:-----|:----------|:------------|
| ComputeAlleleCountsAndHistograms | `data/DataPrep/build/bin` | For each variant it computes NS, AN, AC, AF, Hom, Het, DP, AVGDP, AVGDP_R, AVGGQ, AVGGQ_R, DP_HIST, DP_HIST_R, GQ_HIST, and GQ_HIST_R. Monomorphic variants are dropped (monomorphic variants may arise after subsetting individuals). |
| ComputePercentiles | `data/DataPrep/build/bin` | This program computes percentiles for the QUAL field or any arbitrary numeric INFO field across all variants |
| prepare_sequences.py | `data` | Generates CRAM file with sequences from heterozygous/homozygous samples |

Specify `--help` to see all possible options for each of the above tools (e.g. `Count Alleles --help` or `python prepare_sequences.py --help`).

In addition to the in-house tools listed above, you will need Variant Effect Predictor (VEP) and bcftools.

### Prepare VCF

If applicable, prepare a list of sample ID's that you want to include into BRAVO.
Then, on your raw VCF file(s) run the following:

1. `ComputeAlleleCountsAndHistograms` on your input VCF with raw genotypes:
    `ComputeAlleleCountsAndHistograms -i [input bcf/vcf] -s [samples file] -r [chromosomal region] -o [output vcf.gz] `
2. Variant Effect Predictor (VEP) on the output VCF from step (1).

### Prepare percentiles

To compute percentiles for each INFO field in each variant, use `ComputePercentiles` tool.

### Prepare coverage

1. Use the code in `data/base_coverage/glf2depth/` to create a full coverage file (ie, with coverage for every available base).
Make one `*.full.json.gz` for each chromosome in some directory.

2. Use the scripts in `data/base_coverage/` to bin the coverage.
Make a couple directories with different levels of binning (and again, one `.json.gz` file for each chromosome).

3. Tabix them all.

4. Reference all of the coverage files in `BASE_COVERAGE` in `default.py`.

### Prepare CRAM

BRAVO uses IGV.js to visualize raw sequenced from up to 10 random alternate allele carriers.    

## Load Data

After all data has been prepared and stored in the proper directories, you will need to load that data into the MongoDB database. 

In order to load a human genome into the database, you can either run the `data-import.sh` script directly or use the commands located inside. In either case you will need to modify the script to reflect your environment and needs. By default the script will execute the `INSTALL.sh` script with the default hg38 genome and use 12 threads.

**Warning:** Data initialization may take several hours or longer.

## Data Backup and Restore

Data backup and restoration is handled through MongoDB database dumps to archives. Optionally these can be moved to Google Cloud Storage or the platform of your choice.
Scripts to perform both the database dump as well as the restore can be found in `backup-cmds.sh`. You can modify these scripts to meet the needs of your environment.

**Warning:** Data backup and restoration can take several hours or longer depending on the size of your database, network connection, computational resources, etc..
