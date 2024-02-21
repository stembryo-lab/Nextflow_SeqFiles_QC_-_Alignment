process STAR_INDEX_GENERATE {
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