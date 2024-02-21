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
        --soloCBwhitelist $whitelist \\
        --clipAdapterType CellRanger4 \\
        --soloMultiMappers EM \\
        --soloCellFilter EmptyDrops_CR \\
        --outFilterScoreMin 30 \\
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \\
        --soloUMIfiltering MultiGeneUMI_CR \\
        --soloUMIdedup 1MM_CR \\
        --soloStrand Reverse \\
        --outSAMattributes NH HI AS nM GX GN sM sQ sS '''
    """
}