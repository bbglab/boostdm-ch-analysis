#!/bin/bash

#Prepare access
unset DX_WORKSPACE_ID
dx cd $DX_PROJECT_CONTEXT_ID:

mkdir refs && cd refs/
dx download -r /analysis/refs/vep/
dx download /analysis/refs/hg38_cds.canonical.overlap.gz
dx download /analysis/refs/U2AF1_v*

#Prepare reference samtools
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai 
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.ref_cache.tar.gz
tar xzf Homo_sapiens_assembly38.ref_cache.tar.gz
export REF_PATH="$(pwd)/ref/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s"
export REF_CACHE="$(pwd)/ref/cache/%2s/%2s/%s"

cd vep/
tar xzf homo_sapiens_vep_101_GRCh38.tar.gz

#Prepare environment
cd ../.. && mkdir src && cd src/
dx download -r /analysis/src/U2AF1/*
dx download /analysis/src/cram_UKbio.txt
dx download -r /analysis/src/cram/

cd ../
wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.10.3-Linux-x86_64.sh && bash Miniconda3-py38_4.10.3-Linux-x86_64.sh -b -p /home/dnanexus/miniconda3/
rm Miniconda3-py38_4.10.3-Linux-x86_64.sh
cd /home/dnanexus/miniconda3/bin/
./conda init bash
. /home/dnanexus/.bashrc

cd ~/src/
conda env create -f u2.yml
conda activate u2
