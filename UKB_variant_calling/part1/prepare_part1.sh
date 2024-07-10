#!/bin/bash

#Prepare access
unset DX_WORKSPACE_ID
dx cd $DX_PROJECT_CONTEXT_ID:

#Download references
mkdir refs && cd refs/
dx download /analysis/refs/*

#Prepare reference samtools
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai 
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.ref_cache.tar.gz
tar xzf Homo_sapiens_assembly38.ref_cache.tar.gz
export REF_PATH="$(pwd)/ref/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s"
export REF_CACHE="$(pwd)/ref/cache/%2s/%2s/%s"

#Prepare environment
cd .. && mkdir src && cd src/
dx download /analysis/minimutect_12genes.nf
dx download /analysis/envUKbio.yml
dx download /analysis/nextflow.config
dx download /analysis/bbglab_ukbio.conf
dx download /analysis/run.sh
dx download /analysis/cram_UKbio.txt
dx download -r /analysis/cram/
wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.10.3-Linux-x86_64.sh && bash Miniconda3-py38_4.10.3-Linux-x86_64.sh -b -p /home/dnanexus/miniconda3/
rm Miniconda3-py38_4.10.3-Linux-x86_64.sh
cd /home/dnanexus/miniconda3/bin/
./conda init bash
. /home/dnanexus/.bashrc
cd ~/src/
conda env create -f envUKbio.yml
conda activate nextflow

#Prepare ref
cd ~/refs/
gatk CreateSequenceDictionary -R Homo_sapiens_assembly38.fasta
