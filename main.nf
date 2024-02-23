#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* 
NEXTFLOW PREPROCESSING PIPELINE
-> From .fastq files to .h5ad

## Workflow steps ##
- FASTQC
- TRIMMING
- CREATING GENOME INDEX
- ALIGNMENT
- MTX to H5AD
*/

/* Setting parameters (params.) */
params.timepoint="120h"
params.dataset = "susanne"
params.biosample = "GSM4040776"

/* paths to reads */
params.dataDir = "$baseDir/data/fastqs/${params.timepoint}/${params.dataset}/${params.biosample}"
params.fastqs = "*_{1,2}.fastq.gz"
params.fastqs_trimmed = "*.fq.gz"

params.path2fastqs ="${params.dataDir}/${params.fastqs}"
params.path2fastqs_trimmed ="${params.dataDir}/${params.fastqs_trimmed}"

/* paths to references */
params.annot = "annotation/gencode.vM25.annotation.gtf"
params.fasta = "fasta/GRCm38.p6.genome.fa"
params.whitelist= "whitelist/3M-february-2018.txt"

params.path2ref="$baseDir/references/"
params.outdir="$baseDir/results"

log.info """\
         S C - R N A S E Q - N F   P I P E L I N E    
         ===================================
         # DATASET #
         Timepoint    : ${params.timepoint}
         Dataset      : ${params.dataset}
         Biosample    : ${params.biosample}
                  
         # REFERENCE #
         Genome       : ${params.fasta}
         Annotation   : ${params.annot}

         ===================================
         """
         .stripIndent()

include { FASTQC }                      from './processes/fastqc'
include { MULTIQC }                     from './processes/multiqc'
include { TRIMMING }                    from './processes/trimming'
include { STAR_INDEX_GENERATE }         from './processes/star_index_generate'
include { STAR_ALIGN }                  from './processes/starsolo_align'

workflow {
    'Raw reads channel and QC'
    Channel
        .fromFilePairs(params.path2fastqs, checkIfExists:true)
        .set{ raw_reads_ch }
    fastqc_ch=FASTQC(raw_reads_ch)

    'Trimming using TrimGalore'
    TRIMMING(raw_reads_ch)

    Channel
        .fromPath(params.path2fastqs_trimmed, checkIfExists:true)
        .set { trimmed_reads_ch }

    fastqc_ch=FASTQC(trimmed_reads_ch)
    MULTIQC(fastqc_ch.collect())
    '''
    'Star Genome index generation'
    index_ch=INDEX_GENERATE( params.fasta, 
                             params.annotation )

    STAR_ALIGN( trimmed_reads_ch, 
                index_ch, 
                params.annotation,
                params.whitelist )
'''
}


