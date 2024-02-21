process FASTQC {
    tag "FastQC on $sample_id"
    publishDir "$params.outdir", mode:'copy'
    conda "bioconda::fastqc=0.11.9"

    input:
    tuple val(sample_id), path(reads) 

    output:
    path "${params.dataset}/fastqc/raw/${sample_id}" 

    script:
    """
    mkdir -p "${params.dataset}/fastqc/raw/${sample_id}"

    fastqc -o "${params.dataset}/fastqc/raw/${sample_id}" -q ${reads} 
    """  
}