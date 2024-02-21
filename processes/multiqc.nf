process MULTIQC {
    tag "MultiQC on ${params.dataset} dataset"
    publishDir "$params.outdir", mode:'copy'
    conda "bioconda::multiqc=1.17"

    input:
    path '*'

    output:
    path "${params.dataset}/multiqc/raw"

    script:
    """
    mkdir -p "${params.dataset}/multiqc/raw"

    multiqc . -o "${params.dataset}/multiqc/raw" -n ${params.dataset}_multiqc_report -f -s
    """
}