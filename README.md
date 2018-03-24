Installation
============

1. [System Dependencies](#system-dependencies)
2. [Installation](#installation)
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

Install MongoDB.  Configure it and make it run all the time.

Install python packages in a virtualenv.

    virtualenv venv
    source venv/bin/activate
    pip2 install -r requirements.txt

Some packages will require Python headers (python-dev on some systems).

You probably want to run Flask behind something else, like Apache2.

## Installation

1. Modify the `MONGO` variable in `flask_config.py`, e.g.:
```python
MONGO = {
    'host': 'localhost', # MongoDB hostname
    'port': 27017,       # MongoDB port
    'name': 'bravo'      # Bravo database name
}
```

2. Run `INSTALL.pl` script in `deploy` directory. Specify a human genome build version after `-b`, and a number of threads after `-t`:
```
cd deploy
./INSTALL.pl -b hg38 -t 12
```
Installation may take several hours or longer.

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
    
## Access Control
    
### Authentication

Bravo supports user authentication using Google's OAuth 2.0 protocol, which is optional and is disabled by default. This section describes how to enable it.

First, make sure that your Bravo instance is served using HTTPS protocol.

Second, you need to set up a OAuth with Google. Go [here](https://console.developers.google.com/apis/credentials) and create a project. Your project will get a `Client ID` and a `Client secret`. In the list "Authorized redirect URIs" add your OAuth callback URL, which should look like `https://[bravo base URL]/callback/google` (e.g. `https://mybravo.myinstitution.org/callback/google`).

**Attention!** Don't expose to anyone your `Client ID` and `Client secret`, and make sure you are using HTTPS for your callback URL.

Third, follow these steps to enable authentication in Bravo:
1. Set the `GOOGLE_AUTH` variable in Bravo configuration file to `True`.
2. Assign the `GOOGLE_LOGIN_CLIENT_ID` variable in Bravo configuration file to your `Client ID` from Google.
3. Assign the `GOOGLE_LOGIN_CLIENT_SECRET` variable in Bravo configuration file to your `Client secret` from Google.

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
