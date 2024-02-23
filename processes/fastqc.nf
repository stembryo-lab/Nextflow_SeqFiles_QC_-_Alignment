/* Set output directory */
params.out_fastqc = 'fastqc/'

process FASTQC {
    tag "FastQC on $read_id"
    publishDir "${params.dataDir}/${params.out_fastqc}", mode:'copy'
    conda "bioconda::fastqc=0.11.9"

    input:
    tuple val(read_id), path(reads) 

    output:

    path "*"

    script:

    """

    # create output folder
    mkdir -p ${params.dataDir}/${params.out_fastqc}

    # run FastQC
    fastqc \\
        ${reads} \\
        --outdir ${params.dataDir}/${params.out_fastqc} \\
        --threads $task.cpus

    """  
}

process FASTQC_TRIMMED {
    tag "FastQC on $read_id"
    publishDir "${params.dataDir}/${params.out_fastqc}", mode:'copy'
    conda "bioconda::fastqc=0.11.9"

    input:
    tuple val(read_id), path(reads) 

    output:

    path "*"

    script:

    """

    # create output folder
    mkdir -p ${params.dataDir}/${params.out_fastqc}

    # run FastQC
    fastqc \\
        ${reads} \\
        --outdir ${params.dataDir}/${params.out_fastqc} \\
        --threads $task.cpus

    """  
}