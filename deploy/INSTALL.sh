#!/bin/bash

set -e

RED="\033[0;31m"
GREEN="\033[32m"
NOCOLOR="\033[0m"


download_dbSNP() {
    genome_build=$1
    if [ "$genome_build" == "hg38" ]; then
        echo -e "${GREEN}=> Downloading dbSNP 150 for human genome version ${genome_build}${NOCOLOR}"
        wget -O test.bcp.gz ftp://ftp.ncbi.nlm.nih.gov/snp//organisms/human_9606_b150_GRCh38p7/database/data/organism_data/b150_SNPChrPosOnRef_108.bcp.gz
    elif [ "$genome_build" == "hg19" ]; then 
        echo -e "${GREEN}=> Downloading dbSNP 150 for human genome version ${genome_build}${NOCOLOR}"
        wget -O test.bcp.gz ftp://ftp.ncbi.nlm.nih.gov/snp//organisms/human_9606_b150_GRCh37p13/database/data/organism_data/b150_SNPChrPosOnRef_105.bcp.gz
    else
        echo "Human genome build ${genome_build} is not supported!"
        exit 1
    fi
}

download_GENCODE() {
    genome_build=$1
    if [ "$genome_build" == "hg38" ]; then
        echo -e "${GREEN}=> Downloading GENCODE 27 for human genome version ${genome_build}${NOCOLOR}"
        wget -O test.gtf.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz
    elif [ "$genome_build" == "hg19" ]; then
        echo -e "${GREEN}=> Downloading GENCODE 27 for human genome version ${genome_build}${NOCOLOR}"
        wget -O test.gtf.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/GRCh37_mapping/gencode.v27lift37.annotation.gtf.gz
    else
        echo "Human genome build ${genome_build} is not supported!"
        exit 1
    fi
}

download_dbNSFP() {
    echo -e "${GREEN}=> Downloading dbNSFP 3.5a${NOCOLOR}"
    wget -O dbnsfp.zip ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv3.5a.zip
}

download_canonical_transcripts() {

}

download_OMIM() {

}

echo -e "${RED}DOWNLOADING EXTERNAL DATA${NOCOLOR}"

#download_dbSNP hg19
#download_GENCODE hg19
#download_dbNSFP
