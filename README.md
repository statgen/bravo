Installation
============

1. [System Dependencies](#system-dependencies)
2. [Configuration](#configuration)
3. [Data Preparation](#data-preparation)
   1. [Prepare VCF](#prepare-vcf)
   2. [Prepare percentiles](#prepare-percentiles)
   3. [Prepare coverage](#prepare-coverage)
   4. [Prepare CRAM](#prepare-cram)
4. [Access Control](#access-control)
   1. [Authentication](#authentication)
   2. [Email Whitelist](#email-whitelist)
   3. [Terms Of Use](#terms-of-use)
7. [Google Analytics](#google-analytics)
8. [Start the server](#start-the-server)

## System Dependencies

System dependencies can be installed using the `setup.sh` script.
This script will also create the correct directory structure on the local host.

In addition to these dependencies, it is also recommended to proxy the application behind Apache, Nginx, or similar software.

## Configuration

The configuration for environmental dependent settings can be found at `config/default.py`.
Most settings can be left at their default values, however some will need to be changed to match your environment.


In order to load a human genome into the database, you can either run the `data-import.sh` script directly or use the commands located inside.
In either case you will need to modify the script to reflect your environment and needs. By default the script will execute the `INSTALL.sh` script
with the default hg38 genome and use 12 threads.

**Warning**: Data initialization may take several hours or longer.

## Data Preparation

In the `data/` directory you will find tools/scripts to prepare your data for importing into Mongo database and using in Bravo browser.
Some of these tools are implemented in C++ and need to be compiled as follows:

    cd data/DataPrep/
    cget install .

After successful compilation, the executables will be installed in `data/DataPrep/cget/bin`.
Table below lists all tools that are needed for data preparation.

| Tool | Localtion | Description |
|:-----|:----------|:------------|
| ComputeAlleleCountsAndHistograms | `data/DataPrep/build/bin` | For each variant it computes NS, AN, AC, AF, Hom, Het, DP, AVGDP, AVGGQ, histograms for DP and GQ. Monomorphic variants are dropped (monomorphic variants may arise after subsetting individuals). |
| ComputePercentiles | `data/DataPrep/build/bin` | This program computes percentiles for the QUAL field or any arbitrary numeric INFO field across all variants |
| prepare_sequences.py | `data` | Generates CRAM file with sequences from heterozygous/homozygous samples |

Specify `--help` to see all possible options for each of the above tools (e.g. `Count Alleles --help` or `python prepare_sequences.py --help`).

In addition to the in-house tools listed above, you will need Variant Effect Predictor (VEP) and bcftools.

### Prepare VCF

If applicable, prepare a list of sample ID's that you want to include into Bravo.
Then, on your raw VCF file(s) run the following:

1. `ComputeAlleleCountsAndHistograms` on your input VCF with raw genotypes.
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

Bravo uses IGV.js to visualize raw sequenced from up to 10 random alternate allele carriers.    
    
## Access Control

### Authentication

Bravo supports user authentication using Google's OAuth 2.0 protocol, which is optional and is disabled by default. This section describes how to enable it.

First, make sure that your Bravo instance is served using HTTPS protocol.

Second, you need to set up a OAuth with Google. Go [here](https://console.developers.google.com/apis/credentials) and create a project. Your project will get a `Client ID` and a `Client secret`. In the list "Authorized redirect URIs" add your OAuth callback URL, which should look like `https://[bravo base URL]/callback/google` (e.g. `https://mybravo.myinstitution.org/callback/google`).

**Attention!** Don't expose to anyone your `Client ID` and `Client secret`, and make sure you are using HTTPS for your callback URL.

Third, follow these steps to enable authentication in Bravo:
1. Set the `GOOGLE_AUTH` variable in the Bravo configuration file `config/default.py` to `True`.
2. Assign the `GOOGLE_LOGIN_CLIENT_ID` variable in the Bravo configuration file `config/default.py` to your `Client ID` from Google.
3. Assign the `GOOGLE_LOGIN_CLIENT_SECRET` variable in the Bravo configuration file `config/default.py` to your `Client secret` from Google.

### Email Whitelist
Bravo allows whitelist specific users based on their email address. To enable whitelisting, follow these steps:
1. Set up user authentication as described in [Authentication](#authentication) section.
2. Set the `EMAIL_WHITELIST` variable in Bravo configuration file to `True`.
3. Import list of emails from a text file (one email per line) to the Mongo database:
    
       ./manage.py whitelist -w emails.txt
      
### Terms of Use
If your Bravo users must to agree to any terms/conditions before browsing your data, you need to enable **Terms of Use** page as follows:
1. Set up user authentication as described in [Authentication](#authentication) section.
2. Set the `TERMS` variable in Bravo configuration file to `True`.
3. Write your terms/conditions to the `templates/terms.html` file.


## Google Analytics
This step is optional. Go [here](https://analytics.google.com/analytics/web) and do whatever you have to to get your own `UA-xxxxxx-xx` tracking id.  Put that in `default.py`.  Or just leave the default `UA-01234567-89`, and you won't receive any of the tracking data.

## Start the server

In order to start the service, run `docker-compose up -d` from the projects home directory.

If you are making changed to configuration files or code and need to reload the changes, run `docker-compose up --build -d` in order to rebuild the containers with your updates.
