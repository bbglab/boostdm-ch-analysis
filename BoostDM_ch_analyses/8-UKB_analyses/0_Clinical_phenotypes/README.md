# DOWNLOAD OF UKBIOBANK DATABASE

This folder contains the notebooks to generate the clinical variables in the UKB.
First, clinical data from UKB must be downloaded.

## How to access the data
Acces to UKB sequencing data to carry out the variant calling of blood somatic mutations must be requested via Access Management System (AMS).  
The procedure and conditions to access these datasets are detailed in the following site: https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access

## Steps 

1. Download the "Main dataset" from AMS using the MD5 checksum password

2. Validating the download with "ukbmd5"
```bash
ukbmd5 ukb670124.enc
```

3. Decrypting and uncompressing the dataset with "ukbunpack" 
```bash
ukbunpack ukb670124.enc key_example.key
```
4. Conversion of the dataset: Creation of a dictionary (HTML document)
```bash
ukbconv ukb670124.enc_ukb docs
```

5. Get specific fields related 

    1st. Create file with the data-fields to extract. One data-fields per row.

    2nd. Make conversion using "-i" and "-o":
```bash
ukbconv ukb670124.enc_ukb txt -idata_fields_file -odata_field_output
```