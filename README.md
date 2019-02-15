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

### Prepare VCF

We recommend to process each chromosome separately in parallel. You can further parallelize the process by specifying chromosomal regions in steps (2) and (3).

1. Compile data preparation tools:
   ```
   cd data/DataPrep/
   cget install .
   ```
   After successful compilation, the executables will be installed in `data/DataPrep/cget/bin`.
   
2. Preprare VCF with the INFO fields NS, AN, AC, AF, Hom, Het, DP, AVGDP, AVGDP_R, AVGGQ, AVGGQ_R, DP_HIST, DP_HIST_R, GQ_HIST, and GQ_HIST_R:
   ```
   ./cget/bin/ComputeAlleleCountsAndHistograms -i [input bcf/vcf] -s [samples file] -r [CHR:START-END] -o [output.vcf.gz]
   ```
   Input BCF/VCF must have DP and GQ FORMAT fields. Input BCF/VCF can be accessed both from local and from Google bucket storages. The input samples file (one sample ID per line) and chromosomal region CHR:START-END are optional.

3. Run [Variant Effect Predictor (VEP)](https://www.ensembl.org/vep) on the VCF created in step (2):
   ```
   ./vep -i [input vcf.gz] --plugin LoF[,options]  --assembly [GRCh37/GRCh38] --cache --offline --vcf --sift b --polyphen b --ccds --uniprot --hgvs --symbol --numbers --domains --regulatory --canonical --protein --biotype --af --af_1kg --pubmed --shift_hgvs 0 --allele_number --format vcf --force --buffer_size 100000 --compress_output gzip --no_stats -o [output vcf.gz]
   ```
   Specify [LoF plugin](https://github.com/konradjk/loftee) configuration options as you need.

4. (Optional) Obtain CADD scores from [https://cadd.gs.washington.edu](https://cadd.gs.washington.edu) and annotate VCF from step (3):
   ```
   python add_cadd_scores.py -i [input vcf.gz] -c [cadd_file1.tsv.gz] [cadd_file2.tsv.gz] ...  -o [output vcf.gz]
   ```
   CADD score files must by accompanied by the corresponding index files. If multiple CADD score files are specified, then the maximal CADD score across all files will be used.
<!-- 5. Now you are ready to import VCF's from step (4) into Mongo database. Index all input VCF files with `tabix` and run the following command:
   ```
   python manage.py variants -t [threads] -v [input chr1 vcf.gz] [input chr2 vcf.gz] ...
   ``` -->

### Prepare percentiles
Percentiles must be computed separately for each INFO field.

1. For each VCF INFO field (AVGDP, BQZ, CYZ, DP, FIBC_I, FIBC_P, HWE_SLP_I, HWE_SLP_P, IOR, NM0, NM1, NMZ, QUAL, STZ, SVM, ABE, ABZ) run:
   ```
   ./cget/bin/ComputePercentiles -i [input vcf.gz] -m [INFO field] -t [threads] -f [min MAF] -F [max MAF] -a [allele count] -p [number of perceniles] -d [description] -o [prefix for output files]
   ```
   Examples:
   ```
   ./cget/bin/ComputePercentiles -i /mymachine/myhome/mydata/chr*.mystudy.vcf.gz -t 10 -p 10 -o QUAL
   ./cget/bin/ComputePercentiles -i /mymachine/myhome/mydata/chr*.mystudy.vcf.gz -m ABE -t 10 -p 10 -d "Expected allele Balance towards Reference Allele on Heterozygous Sites" -o ABE
   ```

2. For each INFO field `X` in step (1), you will have two files `X.all_percentiles.json.gz` and `X.variant_percentile.vcf.gz`. The first is a compressed text file with INFO field description and percentiles in JSON format. The second is a compressed VCF file with `X_PCTL` INFO field which stores the corresponding percentile for every variant.
   
3. Index `X.variant_percentile.vcf.gz` using `tabix` and annotate your VCF files from previous step:
   ```
   find . -maxdepth 1 -name "*.variant_percentile..gz" -exec tabix {} \;
   python add_percentiles.py -i [input vcf.gz] - p QUAL.variant_percentile.vcf.gz ABE.variant_percentile.vcf.gz ... -o [output vcf.gz]
   ```

<!-- 3. Import `ALL.all_percentiles.gz` from step (2) into Mongo database:
    ```
    python manage.py metrics -m ALL.all_percentiles.gz
    ```
 4. Update Mongo database with variant percentiles from `*.variant_percentiles.gz` files:
    ```
    [will be added soon]
    ``` -->

### Prepare coverage

To prepare a coverage data for each base-pair position, you can use all your BAM/CRAM files or only a random subset of them (e.g. 1,000) if you need to reduce computational time.

1. For each chromosome and for each BAM/CRAM file extract depth per base-pair:
   ```
   samtools view -q 20 -F 0x0704 -uh [CRAM/BAM file] [chromosome] | samtools calmd -uAEr - [reference FASTA] | bam clipOverlap --in -.ubam --out -.ubam | samtools mpileup -f [reference FASTA] -Q 20 -t DP - | cut -f1-4 | bgzip > [chromosome].[sample].depth.gz
   ```
   In this step we use [clipOverlap from BamUtil](https://genome.sph.umich.edu/wiki/BamUtil:_clipOverlap).
   
2. For each chromosome, create tabix index files for `[chromosome].[sample].depth.gz` e.g.:
   ```
   for f in 10.*.depth.gz; do tabix $f; done
   ```
 
3. For each chromosome, aggregate base-pair coverage acrosss output files `[chromosome].[sample].depth.gz` from step (1):
   ```
   python base_coverage/create_coverage.py -i [files list] aggregate -c [chromosome] -s [start bp] -e [end bp] | bgzip -c > [chromosome].[start].[end].json.gz
   ```
   The `files list` is a text file which lists all output files for a single chromosome from step (1). For example, you can create such list for chromosome 10 with `find . -name "10.*.depth.gz" -printf "%f\n" > depth_files.txt`.
   
   To generate chunks for each chromosome you can use the following command:
   ```
   python base_coverage/create_coverage.py -i [files list] chunk -c [chromosome] -s [chunk size in bp]
   ```
   A typical chunk size is 250,000 bp or 500,000 bp.
   
4. For each chromosome, merge files `[chromosome].[start].[end].json.bgz` from step (3):
   ```
   python base_coverage/merge_coverage.py -i [files list] -o [chromosome].full.json.gz
   ```
   The `files list` is a text file which lists all output files for a single chromosome from step (3).
 
5. After step (4), you should have coverage summary across your samples for each base pair in files `1.full.json.gz`, `2.full.json.gz`, ..., `22.full.json.gz`. For faster web-based visualization, you should prepare several pruned version of the coverage summary e.g.:
   ```
   python base_coverage/prune_coverage.py -i 22.full.json.gz -l 0.25 -o 22.bin_0.25.json.gz
   python base_coverage/prune_coverage.py -i 22.full.json.gz -l 0.50 -o 22.bin_0.50.json.gz
   python base_coverage/prune_coverage.py -i 22.full.json.gz -l 0.75 -o 22.bin_0.75.json.gz
   python base_coverage/prune_coverage.py -i 22.full.json.gz -l 1.00 -o 22.bin_1.00.json.gz
   ```
6. Tabix all coverage summary files.
7. Reference all of the coverage files in `BASE_COVERAGE` in `default.py`.

### Prepare CRAM

BRAVO uses IGV.js to visualize raw sequenced from up to 10 random alternate allele carriers. To enable this visualization BRAVO uses a pre-computed combined CRAM file with all reads from the carriers. We recommend to prepare a separate combined CRAM for each chromosome. The following steps describe hot to prepare the combined CRAM file:

1. In this step, for each variant you will extract IDs for 5 random heterozygous and 5 random homozygous alternate allele carriers from your original VCF file with genotypes. (optional) To speed-up the process, split each chromosome into chunks.
   ```
   RandomHetHom -k 5 -e 1987 --i [input vcf.gz] -s [samples file] -r [CHR:START-END] -o [output carriers_ids.vcf.gz]
   ```
   The value `1987` is a random seed, which you may want to change.
   
2. Prepare a text file (e.g. `samples.txt`) which lists BAM/CRAM file path for each sample:
   ```
   SAMPLEID1    /drive1/batch1/sampleid1.bam
   SAMPLEID2    /drive1/batch1/sampleid2.bam
   SAMPLEID3    /drive1/batch2/sampleid3.bam
   ...
   ```
3. In this step, you will create a combined CRAM file, by extracting all carriers' reads that overlap +/-100bp window around each variant.
   ```
   python prepare_sequences.py cram -i [carriers_ids.vcf.gz] -c samples.txt -w 100 -o [output combined.cram]
   ```

Note: sample IDs and read names in the combined CRAM files will be anonymized.

## Load Data

After all data has been prepared and stored in the proper directories, you will need to load that data into the MongoDB database. 

In order to load a human genome into the database, you can either run the `data-import.sh` script directly or use the commands located inside. In either case you will need to modify the script to reflect your environment and needs. By default the script will execute the `INSTALL.sh` script with the default hg38 genome and use 12 threads.

**Warning:** Data initialization may take several hours or longer.

## Data Backup and Restore

Data backup and restoration is handled through MongoDB database dumps to archives. Optionally these can be moved to Google Cloud Storage or the platform of your choice.
Scripts to perform both the database dump as well as the restore can be found in `backup-cmds.sh`. You can modify these scripts to meet the needs of your environment.

**Warning:** Data backup and restoration can take several hours or longer depending on the size of your database, network connection, computational resources, etc..
