process TRIMMING {

    tag "TrimGalore on ${reads[1]}"
    publishDir "$params.dataDir", mode:'copy'
    conda "bioconda::trim-galore=0.6.7"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    path "trimmed/${sample_id}"

    script:

    """
    mkdir -p "trimmed/${sample_id}"

    cp ${reads[0]} "trimmed/${sample_id}"
    trim_galore ${reads[1]} --illumina -q 20 -o "trimmed/${sample_id}"
    """

}