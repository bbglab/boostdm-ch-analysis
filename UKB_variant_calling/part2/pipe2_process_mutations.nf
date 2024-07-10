#!/usr/bin/env nextflow

//nextflow run pipe2_process_mutations.nf -profile bbglab --in '/home/dnanexus/data/vcftmp/' --out 'test_uploadResults' -resume 

nextflow.enable.dsl=2

params.python1 = '/home/dnanexus/src/Nextflow_filter_ch_variants.py'
params.python2 = '/home/dnanexus/src/Nextflow_filter_VEP_variants.py'

ivcf = Channel.fromPath(params.in+'*')

//GETVCF
process GETVCF {
    //publishDir "results/crams", mode: 'copy'
    input:
    path sample_id

    output:
    path "${sample_id}.filtered.vcf.gz"
    shell:
    '''
    file=$(basename -- !{sample_id} .tmp)
    IFS='/'
    read -a strarr <<< !{params.out}
    ofolder=${strarr[3]}
    dx download "/analysis/output/"$ofolder"/"${file}".filtered.vcf.gz"
    '''
}

process CH_FILTER {

    input:
    path sample_id
    path python1

    output:
    path "${sample_id.simpleName}.maf"
    script:
    """
    python ${python1} -i ${sample_id} -o ./
    """
}
process VEP {
    input:
    path(sample_id)
    output:
    path "${sample_id.simpleName}.vep"
    script:
    """
    vep -i ${sample_id} -o ${sample_id.simpleName}.vep --assembly GRCh38 --no_stats --cache --offline --symbol --protein --vcf --canonical --af_gnomad --af_1kg --dir /home/dnanexus/refs/vep --compress_output gzip
    """
}
process VEP_FILTER {
    input:
    path(sample_id)
    path python2
    output:
    path "${sample_id.simpleName}.filt.gz"
    script:
    """
    python ${python2} -i ${sample_id} -o ./
    """
}

process UPLOAD {
    input:
    path(sample_id)
    script:
    """
    md5sum ${sample_id} > ${sample_id}.md5
    dx upload --singlethread -p --destination ${params.out} ${sample_id}
    dx upload --singlethread -p --destination ${params.out} ${sample_id}.md5
    """
}

workflow {
    GETVCF(ivcf)
    CH_FILTER(GETVCF.out, params.python1)
    VEP(CH_FILTER.out)
    VEP_FILTER(VEP.out, params.python2)
    UPLOAD(VEP_FILTER.out)
}
