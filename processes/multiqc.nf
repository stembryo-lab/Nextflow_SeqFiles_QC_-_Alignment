/* Set output directory */
params.out_multiqc = 'multiqc/'

process MULTIQC {
    tag "MultiQC on ${params.biosample}"
    publishDir "${params.dataDir}/${params.out_multiqc}", mode:'copy'
    conda "bioconda::multiqc=1.17"

    input:
    path '*'

    output:
    path "*"

    script:
    """
    mkdir -p "${params.dataset}/${params.out_multiqc}"

    multiqc . \\
    --outdir "${params.dataDir}/${params.out_multiqc}" \\
    -n ${params.biosample}_report -f -s
    """
}