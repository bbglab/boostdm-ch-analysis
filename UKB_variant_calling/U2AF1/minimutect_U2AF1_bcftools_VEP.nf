#!/usr/bin/env nextflow

//nextflow run minimutect_U2AF1_bcftools_VEP.nf -profile bbglab --in ../data/temporal/wnode1/ --out oo/
nextflow.enable.dsl=2

params.bcfpy = '/home/dnanexus/src/filter_bcf_variants.py'
params.filterpy = '/home/dnanexus/src/Nextflow_filter_VEP_variants.py'
params.ref = '/home/dnanexus/refs/Homo_sapiens_assembly38.fasta'
params.refindex = '/home/dnanexus/refs/Homo_sapiens_assembly38.fasta.fai'
params.refdict = '/home/dnanexus/refs/Homo_sapiens_assembly38.dict'
params.bed1 = '/home/dnanexus/refs/U2AF1_v1.bed'
params.bed2 = '/home/dnanexus/refs/U2AF1_v2.bed'

icrams = Channel.fromFilePairs(params.in+'/*.tmp{,.crai}' )

process GETCRAM {
    input:
    tuple val(sample_id), path(bam)
    output:
    tuple val(sample_id), path("${sample_id}.cram{,.crai}")
    shell:
    '''
    file=$(basename -- !{sample_id} .tmp)
    folder=${file:0:2}
    dx download "/Bulk/Exome\\ sequences/Exome\\ OQFE\\ CRAM\\ files/"$folder"/"$file".cram"
    dx download "/Bulk/Exome\\ sequences/Exome\\ OQFE\\ CRAM\\ files/"$folder"/"$file".cram.crai"
    '''
}

process BCF {
    input:
    tuple val(sample_id), path(cram)
    path ref
    path refindex
    path refdict 
    path bed1
    path bed2
    output:
    tuple val(sample_id), path("${sample_id}_*.bcf.txt")
    script:
    """
    bcftools mpileup --skip-indels -R ${bed1} -f ${ref} --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR ${sample_id}.cram > ${sample_id}_1.bcf.txt
    bcftools mpileup --skip-indels -R ${bed2} -f ${ref} --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR ${sample_id}.cram > ${sample_id}_2.bcf.txt
    """
}

process BCF_FILTER {
    input:
    tuple val(sample_id), path(bcf)
    path bcfpy
    output:
    tuple val(sample_id), path("${sample_id}.maf.gz")
    script:
    """
    python ${bcfpy} -i1 ${sample_id}_1.bcf.txt -i2 ${sample_id}_2.bcf.txt -o ./
    """
}

process VEP {
    input:
    tuple val(sample_id), path(maf)
    output:
    tuple val(sample_id), path("${sample_id}.vep.gz")
    shell:
    """
    lines="\$(zcat !{sample_id}.maf.gz | wc -l)"
    if [ \$lines -gt 1 ]; then
        vep -i !{sample_id}.maf.gz -o !{sample_id}.vep.gz --assembly GRCh38 --no_stats --cache --offline --symbol --protein --vcf --canonical --af_gnomad --af_1kg --dir /home/dnanexus/refs/vep --compress_output gzip
    else
    	echo "empty" | gzip > !{sample_id}.vep.gz;
    fi
    """
}

process VEP_FILTER {
    input:
    tuple val(sample_id), path(vep)
    path filterpy
    output:
    tuple val(sample_id), path("${sample_id}.filt.gz")
    shell:
    """
    lines="\$(zcat !{sample_id}.vep.gz | wc -l)"
    if [ \$lines -gt 1 ]; then
    	python !{filterpy} -i !{sample_id}.vep.gz -o ./;
    else
    	touch !{sample_id}.filt.gz;
    fi
    """


}

process UPLOAD {
    input:
    tuple val(sample_id), path(vepFilt)
    script:
    """
    dx upload --singlethread -p --destination ${params.out} ${sample_id}.filt.gz
    """
}

workflow {
    GETCRAM(icrams)
    BCF(GETCRAM.out, params.ref, params.refindex, params.refdict, params.bed1, params.bed2)
    BCF_FILTER(BCF.out, params.bcfpy)
    VEP(BCF_FILTER.out)
    VEP_FILTER(VEP.out, params.filterpy)
    UPLOAD(VEP_FILTER.out)
}
