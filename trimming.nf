process TRIMMING {

    tag "TrimGalore on ${reads[1]}"
    publishDir "${params.dataDir}", mode:'copy'
    conda "bioconda::trim-galore=0.6.7"
    
    input:
    tuple val(read_id), path(reads) 

    output:
    path "*"

    script:

    """
    trim_galore \\
        ${reads[1]} \\
        --illumina \\
        --quality 20 \\

    """

}