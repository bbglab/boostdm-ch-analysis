#!/bin/bash

#Prepare access
unset DX_WORKSPACE_ID
dx cd $DX_PROJECT_CONTEXT_ID:

mkdir refs && cd refs/
dx download -r /analysis/refs/vep/
dx download /analysis/refs/hg38_cds.canonical.overlap.gz

cd vep/
tar xzf homo_sapiens_vep_101_GRCh38.tar.gz

#Prepare environment
cd ../.. && mkdir src && cd src/
dx download -r /analysis/src/part2/*
dx download /analysis/src/nextflow.config
dx download /analysis/src/cram_UKbio.txt
dx download -r /analysis/src/cram/

wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.10.3-Linux-x86_64.sh && bash Miniconda3-py38_4.10.3-Linux-x86_64.sh -b -p /home/dnanexus/miniconda3/
rm Miniconda3-py38_4.10.3-Linux-x86_64.sh
cd /home/dnanexus/miniconda3/bin/
./conda init bash
. /home/dnanexus/.bashrc
cd ~/src/
conda env create -f envUKbio_part2vep.yml
conda activate nextflowp2

