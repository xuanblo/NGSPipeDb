# ----------------------------------------------------------------------------
# step1 bowtie2 build
genome_index_outdir = join(mapping_outdir, "genome_index")
genome_index_path = join(mapping_outdir, "genome_index", genome_index_prefix)
# ----------------------------------------------------------------------------
rule index_genome_by_bowtie2_build:
    message:
        '''
        ------------------------------
        bowtie2-build
        ------------------------------
        '''
    input:
        genomeFasta = config["genomeFasta"],
    output:
        indexok = touch(join(genome_index_outdir, "index.ok")),
    benchmark:
        join(genome_index_outdir, "benchmark.txt")
    log:
        join(genome_index_outdir, "index.log")
    threads:
        5
    shell:
        '''
        bowtie2-build {input.genomeFasta} {genome_index_path} 1>{log} 2>&1;
        '''

# ----------------------------------------------------------------------------
# step2 align with bowtie2
# ----------------------------------------------------------------------------
rule align_genome_by_bowtie2:
    message:
        '''
        bowtie2 mapping
        '''
    input:
        read1 = join(qc_outdir, "{sample}", "{sample}.cleanR1.fq.gz"),
        read2 = join(qc_outdir, "{sample}", "{sample}.cleanR2.fq.gz"),
        indexok = join(genome_index_outdir, "index.ok"),
    output:
        sorted_bam = join(mapping_outdir, "{sample}", "{sample}.sorted.bam")
    wildcard_constraints:
        dataset="\D+"
    benchmark:
        join(mapping_outdir, "{sample}", "benchmark.txt")
    log:
        join(mapping_outdir, "{sample}", "hisat2.log")
    threads:
        5
    shell:
        '''
        bowtie2 -x {genome_index_path} -p {threads} -1 {input.read1} -2 {input.read2} | samtools sort -o {output.sorted_bam} - 1>{log} 2>&1
        '''


