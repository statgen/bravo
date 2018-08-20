#!/bin/bash
set -e

RED="\033[0;31m"
GREEN="\033[32m"
NOCOLOR="\033[0m"

genome_build=hg38
deploy_dir=${PWD}/temp
log_file=${PWD}/INSTALL.log
threads=1

setup_tabix_tools() {
    echo -e "${GREEN}=> Setup bgzip and tabix${NOCOLOR}"
    start_dir=${PWD}
    { cd ${deploy_dir} \
        && wget https://github.com/samtools/htslib/releases/download/1.7/htslib-1.7.tar.bz2 \
        && tar -xjvf htslib-1.7.tar.bz2 \
        && cd htslib-1.7 \
        && ./configure --prefix=${PWD} \
        && make \
        && make install; } >> ${log_file} 2>&1
    if [ ! $? == 0 ]; then
        echo "Error occured! See ${log_file} for more details."
        exit 1
    fi
    tabix=${PWD}/bin/tabix
    bgzip=${PWD}/bin/bgzip
    cd ${start_dir}
}

setup_ensembl_api() {
    echo -e "${GREEN}=> Setup BioPerl and Ensembl API${NOCOLOR}"
    start_dir=${PWD}
    { cd ${deploy_dir} \
        && wget ftp://ftp.ensembl.org/pub/ensembl-api.tar.gz \
        && wget https://cpan.metacpan.org/authors/id/C/CJ/CJFIELDS/BioPerl-1.6.1.tar.gz \
        && tar -zxf ensembl-api.tar.gz \
        && tar -zxf BioPerl-1.6.1.tar.gz; } >> ${log_file} 2>&1
    if [ ! $? == 0 ]; then
        echo "Error occured! See ${log_file} for more details."
        exit 1
    fi
    PERL5LIB=${PERL5LIB}:${PWD}/BioPerl-1.6.1
    PERL5LIB=${PERL5LIB}:${PWD}/ensembl/modules
    PERL5LIB=${PERL5LIB}:${PWD}/ensembl-compara/modules
    PERL5LIB=${PERL5LIB}:${PWD}/ensembl-variation/modules
    PERL5LIB=${PERL5LIB}:${PWD}/ensembl-funcgen/modules
    export PERL5LIB
    cd ${start_dir}
}

download_GENCODE() {
    start_dir=${PWD}
    if [ "$genome_build" == "hg38" ]; then
        echo -e "${GREEN}=> Downloading GENCODE 27 for human genome version ${genome_build}${NOCOLOR}"
        url=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz
    elif [ "$genome_build" == "hg19" ]; then
        echo -e "${GREEN}=> Downloading GENCODE 27 for human genome version ${genome_build}${NOCOLOR}"
        url=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/GRCh37_mapping/gencode.v27lift37.annotation.gtf.gz
    fi
    { cd ${deploy_dir} \
        && wget -O gencode.gtf.gz ${url}; } >> ${log_file} 2>&1
    if [ ! $? == 0 ]; then
        echo "Error occured! See ${log_file} for more details."
        exit 1
    fi
    gencode_file=${PWD}/gencode.gtf.gz
    cd ${start_dir}
}

download_genenames() {
    echo -e "${GREEN}=> Downloading gene names from HGNC${NOCOLOR}"
    start_dir=${PWD}
    { cd ${deploy_dir} \
        && wget -O hgnc.txt ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt \
        && gzip hgnc.txt; } >> ${log_file} 2>&1
    if [ ! $? == 0 ]; then
        echo "Error occured! See ${log_file} for more details."
        exit 1
    fi  
    hgnc_file=${PWD}/hgnc.txt.gz
    cd ${start_dir}
}

download_OMIM() {
    dbversion=`curl -s "http://useast.ensembl.org/biomart/martservice?type=registry" | grep ENSEMBL_MART_ENSEMBL | grep -oP 'displayName="\K[^"]*'`
    echo -e "${GREEN}=> Downloading OMIM descriptions from BioMart ${dbversion}${NOCOLOR}"
    start_dir=${PWD}
    { cd ${deploy_dir} \
        && wget -O omim.txt 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="0" count="" datasetConfigVersion="0.6"><Dataset name="hsapiens_gene_ensembl" interface="default"><Attribute name="ensembl_gene_id"/><Attribute name="ensembl_transcript_id"/><Attribute name="mim_gene_accession"/><Attribute name="mim_gene_description"/></Dataset></Query>' \
        && gzip omim.txt; } >> ${log_file} 2>&1
    if [ ! $? == 0 ]; then
        echo "Error occured! See ${log_file} for more details."
        exit 1
    fi
    omim_file=${PWD}/omim.txt.gz
    cd ${start_dir}
}

download_dbSNP() {
    start_dir=${PWD}
    if [ "$genome_build" == "hg38" ]; then
        echo -e "${GREEN}=> Downloading dbSNP 150 for human genome version ${genome_build}${NOCOLOR}"
        url="ftp://ftp.ncbi.nlm.nih.gov/snp//organisms/human_9606_b150_GRCh38p7/database/data/organism_data/b150_SNPChrPosOnRef_108.bcp.gz"
    elif [ "$genome_build" == "hg19" ]; then 
        echo -e "${GREEN}=> Downloading dbSNP 150 for human genome version ${genome_build}${NOCOLOR}"
        url="ftp://ftp.ncbi.nlm.nih.gov/snp//organisms/human_9606_b150_GRCh37p13/database/data/organism_data/b150_SNPChrPosOnRef_105.bcp.gz"
    fi
    { cd ${deploy_dir} \
        && wget -O dbsnp.bcp.gz ${url} \
        && gzip -dc dbsnp.bcp.gz | awk '{OFS = "\t"; if ($3 != "") { print $1,$2,$3 }}' | sort -k2,2 -k3,3n | ${bgzip} -c > dbsnp.tsv.gz \
        && ${tabix} -s 2 -b 3 -e 3 dbsnp.tsv.gz \
        && rm dbsnp.bcp.gz; } >> ${log_file} 2>&1
    if [ ! $? == 0 ]; then
        echo "Error occured! See ${log_file} for more details."
        exit 1
    fi
    dbsnp_file=${PWD}/dbsnp.tsv.gz
    cd ${start_dir}
}

download_canonical_transcripts() {
    # we should not worry about human genome build version here, because we don't need genomic coordinates
    # we want just a specific version of Ensembl
    dbversion="homo_sapiens_core_91_38"
    echo -e "${GREEN}=> Downloading canonical transcripts from Ensembl 91${NOCOLOR}"
    start_dir=${PWD}
    { cd ${deploy_dir} \
        && perl ${start_dir}/download_canonical_transcripts.pl homo_sapiens_core_91_38 | gzip -c > canonical_transcripts.txt.gz; } >> ${log_file} 2>&1
    if [ ! $? == 0 ]; then
        echo "Error occured! See ${log_file} for more details."
        exit 1
    fi  
    canonical_transcripts_file=${PWD}/canonical_transcripts.txt.gz
    cd ${start_dir}
}


while getopts ":hb:d:t:c:" opt; do
    case $opt in
        h)  
            echo "This script will setup Bravo."
            echo ""
            echo "./INSTALL.sh"
            echo -e "\t-h              Displays this help message."
            echo -e "\t-b build        Human genome build version: hg38 (default), hg19."
            echo -e "\t-d directory    Working directory for temporary deployment files and scripts. Default is ./temp."
            echo -e "\t-t threads      Number of parallel threads. Default is 1."
            echo ""
            exit 0
            ;;
        b)
            genome_build=$OPTARG
            ;;
        d)
            working_dir=$OPTARG
            ;;
        t)
            threads=$OPTARG
            ;;
        \?)
            echo "Invalid option: -$OPTARG"
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires and argument."
            exit 1
            ;;
    esac
done

if [ "$genome_build" != "hg38" ] && [ "$genome_build" != "hg19" ]; then
    echo "Genome build $genome_build is not supported."
    exit 1
fi

if [ -d "${deploy_dir}" ]; then
    echo "Temporary deployment directory ${deploy_dir} already exists. Please remove it or provide an alternative."
    exit 1
fi

mkdir ${deploy_dir}

echo "Setting up Bravo for human genome version ${genome_build}."
echo "Writing temporary deployment files to ${deploy_dir}."
echo "Writing log information to ${log_file}."
echo -n "" > ${log_file}

echo -e "${RED}DOWNLOADING REQUIRED TOOLS/LIBRARIES${NOCOLOR}"
command -v tabix > /dev/null && command -v bgzip > /dev/null
if [ ! $? == 0 ]; then
    setup_tabix_tools
else
    tabix=`command -v tabix`
    bgzip=`command -v bgzip`
fi
setup_ensembl_api

echo -e "${RED}DOWNLOADING EXTERNAL DATA${NOCOLOR}"
download_GENCODE
download_genenames
download_OMIM
download_dbSNP
download_canonical_transcripts

echo -e "${RED}LOADING EXTERNAL DATA TO DATABASE${NOCOLOR}"

python ../manage.py genes -t ${canonical_transcripts_file} -m ${omim_file} -f ${hgnc_file} -g ${gencode_file}
python ../manage.py dbsnp -d ${dbsnp_file} -t ${threads}

rm -rf ${deploy_dir}

echo -e "${RED}DONE${NOCOLOR}"
