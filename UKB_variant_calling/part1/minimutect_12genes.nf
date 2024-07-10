#!/usr/bin/env nextflow

//nextflow run minimutect_12genes.nf -profile bbglab --in '/home/dnanexus/data/tmp5/' --out 'test_uploadResults' -resume 
//renamed from pipe_santid_UKbio_CheckOut_getfiles.nf

nextflow.enable.dsl=2

params.runid = 'run_uploadIncluded'

params.bed = '/home/dnanexus/refs/CHgenes12_coordinates.bed'
params.germline_resource = '/home/dnanexus/refs/somatic-hg38_af-only-gnomad.hg38.vcf.gz'
params.germline_resourceIndex = '/home/dnanexus/refs/somatic-hg38_af-only-gnomad.hg38.vcf.gz.tbi'
params.path_small_exac = '/home/dnanexus/refs/somatic-hg38_small_exac_common_3.hg38_selectedCHROM.vcf.gz'
params.path_small_exacIndex = '/home/dnanexus/refs/somatic-hg38_small_exac_common_3.hg38_selectedCHROM.vcf.gz.tbi'
params.ref = '/home/dnanexus/refs/Homo_sapiens_assembly38.fasta'
params.refindex = '/home/dnanexus/refs/Homo_sapiens_assembly38.fasta.fai'
params.refdict = '/home/dnanexus/refs/Homo_sapiens_assembly38.dict'


icrams = Channel.fromFilePairs(params.in+'/*.tmp{,.crai}' )

process GETCRAM {
    input:
    tuple val(sample_id), path(bam)
    output:
    tuple val(sample_id), path("${sample_id}.cram")
    shell:
    '''
    file=$(basename -- !{sample_id} .tmp)
    folder=${file:0:2}
    dx download "/Bulk/Exome\\ sequences/Exome\\ OQFE\\ CRAM\\ files/"$folder"/"$file".cram"
    dx download "/Bulk/Exome\\ sequences/Exome\\ OQFE\\ CRAM\\ files/"$folder"/"$file".cram.crai"
    '''
}

process MINIBAM {
    input:
    tuple val(sample_id), path(bam)
    path bed 
    output:
    tuple val(sample_id), path("${sample_id}.mini.bam")
    script:
    """
    samtools view --threads 1 -b -L ${bed} ${sample_id}.cram -o ${sample_id}.minitemp.bam
    picard BuildBamIndex I=${sample_id}.minitemp.bam
    samtools view -h ${sample_id}.minitemp.bam | grep -P -v '@SQ\tSN:chr3|@SQ\tSN:chr[5-9]|@SQ\tSN:chr1[0-4]|@SQ\tSN:chr16|@SQ\tSN:chr18|@SQ\tSN:chr19|@SQ\tSN:chrUn|@SQ\tSN:HLA|_alt\tLN|_random\tLN|@SQ\tSN:chr[A-Z]' > ${sample_id}.mini.sam
    samtools view -S -b  ${sample_id}.mini.sam >  ${sample_id}.mini.bam
    samtools index ${sample_id}.mini.bam
    """
}

process MUTECT2 {
    input: 
    tuple val(sample_id), path(minibam)
    path ref
    path refindex
    path germline_resource
    path germline_resourceIndex
    path dict
    path bed
    output:
    tuple val(sample_id), path("${sample_id}*.gz*")
    script:
    """
    gatk Mutect2 --java-options "-Xmx16G" -R ${ref}  -I ${sample_id}.mini.bam  -tumor ${sample_id} --germline-resource ${germline_resource} -O ${sample_id}.vcf.gz --callable-depth 3 --f1r2-tar-gz ${sample_id}_f1r2.tar.gz -L ${bed}
    """
}

process ORIENTATION_MODEL {
    input:
    tuple val(sample_id), path(f1r2)
    output:
    tuple val(sample_id), path("${sample_id}_read-orientation-model.tar.gz")
    script:
    """
    gatk LearnReadOrientationModel -I ${sample_id}_f1r2.tar.gz -O ${sample_id}_read-orientation-model.tar.gz
    """
}

process PILEUP {
    input:
    tuple val(sample_id), path(minibam)
    path small_exac
    path small_exac_index
    output:
    tuple val(sample_id), path("${sample_id}_getpileupsummaries.table")
    script:
    """
    gatk GetPileupSummaries -I ${sample_id}.mini.bam -V ${small_exac} -L ${small_exac} -O ${sample_id}_getpileupsummaries.table
    """
}

process CONTAMINATION {
    input:
    tuple val(sample_id), path(summaries)
    output:
    //tuple val(sample_id), path("${sample_id}_segments_cl.table"), emit: contamination_segments
    //tuple val(sample_id), path("${sample_id}_contamination_cl.table"), emit: contamination_table
    tuple val(sample_id), path("${sample_id}_*cl.table")
    script:
    """
    gatk CalculateContamination  -I ${sample_id}_getpileupsummaries.table -tumor-segmentation ${sample_id}_segments_cl.table -O ${sample_id}_contamination_cl.table
    """
}

process FILTER {
    input:
    tuple val (sample_id), path(vcf_unfilt), path(orientation), path(contaout)
    path ref
    path refindex
    path dict
    output:
    tuple val(sample_id), path("${sample_id}.filtered.vcf.gz")
    script:
    """
    gatk FilterMutectCalls -V ${sample_id}.vcf.gz --tumor-segmentation  ${sample_id}_segments_cl.table -R ${ref} --contamination-table ${sample_id}_contamination_cl.table --ob-priors ${sample_id}_read-orientation-model.tar.gz --max-events-in-region 3 -O ${sample_id}.filtered.vcf.gz
    """
}

process UPLOAD {
    input:
    tuple val (sample_id), path(cram), path(mbam), path(vcf_unfilt), path(orientation), path(contaout), path(vcf_filt)
    script:
    """
    md5sum ${sample_id}.filtered.vcf.gz > ${sample_id}.filtered.vcf.gz.md5
    dx upload --wait --singlethread --no-progress --parents --path '${params.out}' ${sample_id}.filtered.vcf.gz
    dx upload --wait --singlethread --no-progress --parents --path '${params.out}' ${sample_id}.filtered.vcf.gz.md5
    """
}


workflow {
    GETCRAM(icrams)

    MINIBAM(GETCRAM.out, params.bed)

    MUTECT2(MINIBAM.out, params.ref, params.refindex, params.germline_resource, params.germline_resourceIndex, params.refdict, params.bed)

    ORIENTATION_MODEL(MUTECT2.out)

    PILEUP(MINIBAM.out,params.path_small_exac,params.path_small_exacIndex)

    CONTAMINATION(PILEUP.out)

    FILTER(MUTECT2.out.join(ORIENTATION_MODEL.out).join(CONTAMINATION.out),params.ref, params.refindex, params.refdict)

    UPLOAD(GETCRAM.out.join(MINIBAM.out).join(MUTECT2.out).join(ORIENTATION_MODEL.out).join(CONTAMINATION.out).join(FILTER.out))
}
