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


/* 1. Performing QC on fastq files */

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

/* Merged report with MultiQC */

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

/* Trimming adapters and QC on cDNA reads */

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

process FASTQC_TRIMMED {
    tag "FastQC on trimmed $sample_id"
    publishDir "$params.outdir", mode:'copy'
    conda "bioconda::fastqc=0.11.9"

    input:
    tuple val(sample_id), path(reads) 

    output:
    path "${params.dataset}/fastqc/trimmed/${sample_id}/" 

    script:
    """
    mkdir -p "${params.dataset}/fastqc/trimmed/${sample_id}"

    fastqc -o "${params.dataset}/fastqc/trimmed/${sample_id}/" -q ${reads} 
    """  
}

process MULTIQC_TRIMMED {
    tag "MultiQC on trimmed and untrimmed ${params.dataset}"
    publishDir "$params.outdir", mode:'copy'
    conda "bioconda::multiqc=1.17"

    input:
    path '*'

    output:
    path "${params.dataset}/multiqc/trimmed"

    script:
    """
    mkdir -p "${params.dataset}/multiqc/trimmed"

    multiqc . -o "${params.dataset}/multiqc/trimmed" -n ${params.dataset}_trimmed_multiqc_report -f -s
    """
}

process INDEX_GENERATE {
    tag 'Star genome index generation'
    publishDir "./references"
    conda "bioconda::star=2.7.10a bioconda::samtools=1.16.1 conda-forge::gawk=5.1.0"

    input:
    path fasta
    path annot
     
    output:
    path 'genome_index' 

    script:       
    """
    mkdir genome_index
    
    STAR \
    --runMode genomeGenerate \
    --runThreadN 4 \
    --genomeDir genome_index \
    --genomeFastaFiles ${fasta} \
    --sjdbGTFfile ${annot} \
    --genomeSAindexNbases 12 \
    --genomeSAsparseD 3 \
    --limitGenomeGenerateRAM 15000000000

    """
}

'   --sjdbOverhang read_length-1'

process STAR_ALIGN {
    tag "Alignment on $sample_id"
    publishDir "$params.outdir"
    conda "bioconda::star=2.7.10a bioconda::samtools=1.16.1 conda-forge::gawk=5.1.0"
    
    input:
    //
    // Input reads are expected to come as: [ meta, [ pair1_read1, pair1_read2, pair2_read1, pair2_read2 ] ]
    // Input array for a sample is created in the same order reads appear in samplesheet as pairs from replicates are appended to array.
    //
    tuple val(sample_id), path(reads)
    path  index
    path  gtf
    path  whitelist

    output:
    path 'alignment/${params.dataset}/${sample_id}'
    tuple val(sample_id), path('*d.out.bam')      
    tuple val(sample_id), path('*.Solo.out')     
    tuple val(sample_id), path('*Log.final.out')  
    tuple val(sample_id), path('*Log.out')         
    tuple val(sample_id), path('*Log.progress.out')

    script:

    """
    mkdir -p alignment/${params.dataset}/${sample_id}

    STAR \\
        --genomeDir $index \\
        --genomeLoad LoadAndKeep \\
        --readFilesIn ${reads[1]} ${reads[0]} \\
        --runThreadN 8 \\
        --outFileNamePrefix alignment/${params.dataset}/${sample_id}/ \\
        --soloType CB_UMI_Simple \\
        --soloFeatures GeneFull \\
        --soloCBstart 1 \\
        --soloCBlen 16 \\
        --soloUMIstart 17 \\
        --soloUMIlen 12 \\

    """

}

'''       --soloCBwhitelist $whitelist \\
        --clipAdapterType CellRanger4 \\
        --soloMultiMappers EM \\
        --soloCellFilter EmptyDrops_CR \\
        --outFilterScoreMin 30 \\
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \\
        --soloUMIfiltering MultiGeneUMI_CR \\
        --soloUMIdedup 1MM_CR \\
        --soloStrand Reverse \\
        --outSAMattributes NH HI AS nM GX GN sM sQ sS '''

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


