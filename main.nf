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

params.dataset = ""

params.dataDir = "$baseDir/data/${params.dataset}"
params.read_cbc = "1"
params.read_cdna = "2"
params.raw_reads = "$baseDir/data/${params.dataset}/raw/*_{${params.read_cbc},${params.read_cdna}}.fastq.gz"
params.trimmed_reads = "$baseDir/data/${params.dataset}/trimmed/SRR*/*_{${params.read_cbc},${params.read_cdna}}*.f*.gz"

params.annotation = "$baseDir/references/annot/gencode.vM25.primary_assembly.annotation.gtf"
params.fasta = "$baseDir/references/fasta/GRCm38.primary_assembly.genome.fa"
params.whitelist= "$baseDir/references/whitelist/3M-february-2018.txt"

params.outdir="$baseDir/results"

nf-core modules install fastqc
nf-core modules install trimgalore

log.info """\
         S C - R N A S E Q - N F   P I P E L I N E    
         ===================================
         dataset      : ${params.dataset}
         genome       : ${params.fasta}
         annotation   : ${params.annotation}
         reads        : ${params.raw_reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()

include { FASTQC }                      from '../processes/fastqc'
include { MULTIQC }                     from '../processes/multiqc'
include { TRIMMING }                    from '../processes/trimming'
include { STAR_INDEX_GENERATE }         from '../processes/star_index_generate'
include { STAR_ALIGN }                  from '../processes/star_align'

// TO BE REMOVED IF POSSIBLE
// process FASTQC_TRIMMED {
//     tag "FastQC on trimmed $sample_id"
//     publishDir "$params.outdir", mode:'copy'
//     conda "bioconda::fastqc=0.11.9"

//     input:
//     tuple val(sample_id), path(reads) 

//     output:
//     path "${params.dataset}/fastqc/trimmed/${sample_id}/" 

//     script:
//     """
//     mkdir -p "${params.dataset}/fastqc/trimmed/${sample_id}"

//     fastqc -o "${params.dataset}/fastqc/trimmed/${sample_id}/" -q ${reads} 
//     """  
// }

// TO BE REMOVED IF POSSIBLE
// process MULTIQC_TRIMMED {
//     tag "MultiQC on trimmed and untrimmed ${params.dataset}"
//     publishDir "$params.outdir", mode:'copy'
//     conda "bioconda::multiqc=1.17"

//     input:
//     path '*'

//     output:
//     path "${params.dataset}/multiqc/trimmed"

//     script:
//     """
//     mkdir -p "${params.dataset}/multiqc/trimmed"

//     multiqc . -o "${params.dataset}/multiqc/trimmed" -n ${params.dataset}_trimmed_multiqc_report -f -s
//     """
// }

workflow {
    'Raw reads channel and QC'
    Channel
        .fromFilePairs(params.raw_reads, checkIfExists:true)
        .set{ raw_reads_ch }
    fastqc_ch=FASTQC(raw_reads_ch)
    MULTIQC(fastqc_ch.collect())

    'Trimming using TrimGalore'
    TRIMMING(raw_reads_ch)

    'Trimmed reads channel and QC'
    Channel
        .fromFilePairs(params.trimmed_reads, checkIfExists:true)
        .set { trimmed_reads_ch }
    fastqc_all_ch=FASTQC_TRIMMED(trimmed_reads_ch)
    MULTIQC_TRIMMED(fastqc_all_ch.collect())

    'Star Genome index generation'
    index_ch=INDEX_GENERATE( params.fasta, 
                             params.annotation )

    STAR_ALIGN( trimmed_reads_ch, 
                index_ch, 
                params.annotation,
                params.whitelist )

}


