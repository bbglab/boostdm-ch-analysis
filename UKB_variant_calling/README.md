# UKB variant calling

This folder contains the code to perform the variant calling in the UKB.
It is important to remark that all this analysis were run on the [DNAnexus](https://www.dnanexus.com/) cloud platform.

## How to access the data
Acces to UKB sequencing data to carry out the variant calling of blood somatic mutations must be requested via Access Management System (AMS).  
The procedure and conditions to access these datasets are detailed in the following site: https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access

# Content
- ```crams```: cram file's name list
- ```part1```: variant calling for the 12 genes with boostDM models
    - **bbglab_ukbio.conf**, configuration file containing CPU and memory per subprocess
    - **cram_UKbio.txt**, cram file's name list
    - **envUKbio.yml**, conda environment
    - **minimutect_12genes.nf**, nextflow to perform variant calling
    - **nextflow.config**, nextflow configurations
    - **prepare_part1.sh**, Spcript to prepare environment in a virtual machine
    - **run.sh**, script to run nextflow script
- ```part2```: anotation and filtering
    - **bbglab_ukbio.conf**, configuration file containing CPU and memory per subprocess
    - **cram_UKbio.txt**, cram file's name bam
    - **envUKbio_part2vep.yml**, conda environment
    - **merge.sh**, script to merge final outputs
    - **nextflow.config**, nextflow configurations
    - **Nextflow_filter_ch_variants.py**, python script to filter variants
    - **Nextflow_filter_VEP_variants.py**, python script to filter variants afted VEP annotation
    - **pipe2_process_mutations.nf**, nextflow to perform annotationa and filtering
    - **prepare_part2.sh**, Spcript to prepare environment in a virtual machine
    - **run.sh**, script to run nextflow script
- ```U2AF1```: hause-made pileup to call U2AF1 variants
    - **bbglab_ukbio.conf**, configuration file containing CPU and memory per subprocess
    - **cram_UKbio.txt**, cram file's name bam
    - **run_bcftools.sh**, script to merge final outputs
    - **nextflow.config**, nextflow configurations
    - **filter_bcf_variants.py**, python script to filter variants
    - **Nextflow_filter_VEP_variants.py**, python script to filter variants afted VEP annotation
    - **minimutect_U2AF1_bcftools_VEP.nf**, nextflow to perform annotationa and filtering
    - **prepare_U2AF1.sh**, Spcript to prepare environment in a virtual machine
    - **run.sh**, script to run nextflow script
    - **u2.yml**, singularity image to run VEP
- ```refs```: references necessary for the variant calling

# Dependencies (refs folder)
In the refs folder there are the following files:

1. Bed files
    - **CHgenes12_coordinates.bed'**, bed file of the coordinates of the 12 genes
    - **U2AF1_v1.bed'**, bed file with the coordinates of U2AF1 ENST00000291552 transcript
    - **U2AF1_v2.bed'**, bed file with the coordinates of U2AF1 ENSG00000275895 transcript

2. Germline information
    - **somatic-hg38_af-only-gnomad.hg38.vcf.gz**, : https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
    - **somatic-hg38_af-only-gnomad.hg38.vcf.gz.tbi**, : https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi
    - **somatic-hg38_small_exac_common_3.hg38_selectedCHROM.vcf.gz**, : modified from https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz
    - **somatic-hg38_small_exac_common_3.hg38_selectedCHROM.vcf.gz.tbi**: modified from https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi

3. Reference genome
    In the prepare*.sh file includes code to download it:
    - **Homo_sapiens_assembly38.fasta**, wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
    - **Homo_sapiens_assembly38.fasta.fai**, wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai 
    - **Homo_sapiens_assembly38.dict**, wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.ref_cache.tar.gz

4. Other resources
    - **hg38_cds.canonical.overlap.gz**, list of canonical transcripts
    - **homo_sapiens_vep_101_GRCh38.tar.gz**, downloaded from ensamble_VEP https://www.ensembl.org/info/docs/tools/vep/index.html
    - **vep.simg**, singularity image to run VEP.101 (needs to be created)


# Pipeline 

To perform the variant calling run first scripts in ```part1``` to run Mutect2, ```part2``` to annotate and filter muatations and finally ```U2AF1``` to call mutations in the U2AF1 gene. For every part follow the following steps:

**1. Prepare VM on DNAnexus to run the variant caller**
Create a VM with very little resources.
Copy the scripts folder to the running instance:

```bash
prepare*.sh
```

**2. Execute the script in the VM**
Run the run.sh script but before modify the output folder and check that all the reference data needed have been downloaded in the VM.

```bash
run.sh
```

