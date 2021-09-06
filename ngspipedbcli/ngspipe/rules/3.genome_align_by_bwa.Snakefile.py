# ----------------------------------------------------------------------------
# step1 bowtie2 build
genome_index_outdir = join(mapping_outdir, "genome_index")
genome_index_path = join(mapping_outdir, "genome_index", genome_index_prefix)
# ----------------------------------------------------------------------------
rule index_genome_by_bwa_index:
    message:
        '''
        ------------------------------
        bwa-build
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
        bwa index -p {genome_index_path} {input.genomeFasta} 1>{log} 2>&1;
        '''

# ----------------------------------------------------------------------------
# step2 align with bowtie2
# ----------------------------------------------------------------------------
rule align_genome_by_bwa:
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
        join(mapping_outdir, "{sample}", "bwa.log")
    threads:
        5
    shell:
        '''
        bwa mem -t {threads} {genome_index_path} {input.read1} {input.read2} 2>{log} | samtools sort -o {output.sorted_bam} - 1>>{log} 2>&1;
        '''


