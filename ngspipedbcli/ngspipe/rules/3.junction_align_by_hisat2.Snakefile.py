# ----------------------------------------------------------------------------
# step1 hisat2 build
genome_index_outdir = join(junction_align_outdir, "genome_index")
genome_index_path = join(junction_align_outdir, "genome_index", genome_index_prefix)
# ----------------------------------------------------------------------------
rule index_genome_by_hisat2_build:
    message:
        '''
        ------------------------------
        hisat2-build
        ------------------------------
        '''
    input:
        genomeFasta = config["genomeFasta"],
        genomeAnno = config["genomeAnno"]
    output:
        exon_info = join(genome_index_outdir, "genome.exon"),
        splice_sites = join(genome_index_outdir, "genome.ss")
    benchmark:
        join(genome_index_outdir, "benchmark.txt")
    log:
        join(genome_index_outdir, "index.log")
    conda:
        "../envs/hisat2.yaml"
    threads:
        5
    shell:
        '''
        hisat2_extract_exons.py {input.genomeAnno} > {output.exon_info} 2>>{log};
        hisat2_extract_splice_sites.py {input.genomeAnno} > {output.splice_sites} 2>>{log};
        hisat2-build -p {threads} --ss {output.splice_sites} --exon {output.exon_info} {input.genomeFasta} {genome_index_path} 1>>{log} 2>&1;
        '''

# ----------------------------------------------------------------------------
# step2 align with hisat2
# ----------------------------------------------------------------------------
rule align_genome_by_hisat2:
    message:
        '''
        hisat2
        '''
    input:
        read1 = join(qc_outdir, "{sample}", "{sample}.cleanR1.fq.gz"),
        read2 = join(qc_outdir, "{sample}", "{sample}.cleanR2.fq.gz"),
        splice_sites = join(genome_index_outdir, "genome.ss")
    output:
        sorted_bam = join(junction_align_outdir, "{sample}", "{sample}.sorted.bam")
    wildcard_constraints:
        dataset="\D+"
    benchmark:
        join(junction_align_outdir, "{sample}", "benchmark.txt")
    log:
        join(junction_align_outdir, "{sample}", "hisat2.log")
    threads:
        5
    shell:
        '''
        hisat2 --dta {rna_library} -p {threads} -x {genome_index_path} -1 {input.read1} -2 {input.read2} 2>{log} | samtools sort -@ {threads} -o {output.sorted_bam} 1>>{log} 2>&1
        '''


