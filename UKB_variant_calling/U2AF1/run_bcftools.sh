#!/bin/bash

#CONFIGURATION
#output at RAP
OFOLDER='/analysis/output/bcftools/wnode1/'
#temporal folder for batches
TFOLDER='/home/dnanexus/data/bcftools/temporal/wnode1'
IFILE='cram/cram_UKbio_01.txt'

echo "Checking samples analysed..."
#Clean input folder and set up counters:
COUNTER=0
mkdir -p $TFOLDER
rm -f $TFOLDER/*
#Check files already in results folder
dx ls "$OFOLDER"*1.bcf.txt | sed 's/[.].*$//' > filesDone.txt
xlines=$(wc -l < filesDone.txt)
ylines=$(wc -l < $IFILE)
echo "Number of samples already found at RAP: "$xlines
echo "Samples pending: "$((ylines-xlines))
#Set total counter
TCOUNTER=xlines
#Run loop
echo "Starting analysis..."
while read line; do if [ $COUNTER == 1000 ]; then nextflow run minimutect_U2AF1_bcftools.nf -profile bbglab --in $TFOLDER/ --out $OFOLDER -resume; COUNTER=0; rm $TFOLDER/*; echo "Total samples processed: "$TCOUNTER; fi; line=${line##*/}; s="${line%.cram}"; if ! grep -Fxq $s filesDone.txt; then touch $TFOLDER/$s.tmp; touch $TFOLDER/$s.tmp.crai; (( COUNTER++ )); (( TCOUNTER++ )); fi; done < $IFILE
#and last iteration
if [ $COUNTER != 0 ]; then nextflow run minimutect_U2AF1_bcftools.nf -profile bbglab --in $TFOLDER/ --out $OFOLDER -resume; COUNTER=0; rm $TFOLDER/*; fi;
