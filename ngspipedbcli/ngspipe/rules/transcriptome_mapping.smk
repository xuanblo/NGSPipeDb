

junction_align_method = config['junction_align_method'] # star
junction_align_outdir = join(config["resultsDir"], "mapping")
genome_index_prefix = "genome"
# rna-seq sequencing type, can be fr-firststrand, none, fr-secondstrand
#rna_library = "" # "--rna-strandness RF"(fr-firststrand) or "--rna-strandness FR"(fr-secondstrand)

rna_library = '' if config['junction_align_rna_library'] == None else config['junction_align_rna_library']

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
        exon_info = join(junction_align_outdir, 'hisat2', 'index', "genome.exon"),
        splice_sites = join(junction_align_outdir, 'hisat2', 'index', "genome.ss")
    benchmark:
        join(junction_align_outdir, 'hisat2', 'index', "benchmark.txt")
    log:
        join(junction_align_outdir, 'hisat2', 'index', "index.log.txt")
    conda:
        "../envs/hisat2.yaml"
    threads:
        5
    shell:
        '''
        hisat2_extract_exons.py {input.genomeAnno} > {output.exon_info} 2>>{log};
        hisat2_extract_splice_sites.py {input.genomeAnno} > {output.splice_sites} 2>>{log};
        hisat2-build -p {threads} --ss {output.splice_sites} --exon {output.exon_info} {input.genomeFasta} {junction_align_outdir}/hisat2/index/{genome_index_prefix} 1>>{log} 2>&1;
        '''

# ----------------------------------------------------------------------------
# step2 align with hisat2
# ----------------------------------------------------------------------------
rule transcriptome_align_by_hisat2:
    message:
        '''
        hisat2
        '''
    input:
        read1 = join(qc_outdir, qc_method, "{sample}", "{sample}.cleanR1.fq.gz"),
        read2 = join(qc_outdir, qc_method,"{sample}", "{sample}.cleanR2.fq.gz"),
        splice_sites = join(junction_align_outdir, 'hisat2', 'index', "genome.ss")
    output:
        sorted_bam = join(junction_align_outdir, 'hisat2', 'align', "{sample}", "{sample}.sorted.bam")
    wildcard_constraints:
        dataset="\D+"
    benchmark:
        join(junction_align_outdir, 'hisat2', 'align', "{sample}", "benchmark.txt")
    log:
        join(junction_align_outdir, 'hisat2', 'align', "{sample}", "hisat2.log.txt")
    threads:
        5
    shell:
        '''
        hisat2 --dta {rna_library} -p {threads} -x {junction_align_outdir}/hisat2/index/{genome_index_prefix} -1 {input.read1} -2 {input.read2} 2>{log} | samtools sort -@ {threads} -o {output.sorted_bam} 1>>{log} 2>&1
        '''

rule genome_index_by_star:
    message:
        "Running STAR Alignment index" 
    input:
        genomeFasta = config["genomeFasta"],
        genomeAnno = config["genomeAnno"]
    output:
        index_ok = touch(join(junction_align_outdir, 'star', 'index', 'index.ok'))
    threads: 8
    log:
        log_file = join(junction_align_outdir, 'star', 'index', "index.log.txt")
    priority: 10
    benchmark: join(junction_align_outdir, 'star', 'index', "benchmark.txt")
    params:
        outdir = join(junction_align_outdir, 'star', 'index')
    shell:
        '''
        STAR --runMode genomeGenerate --runThreadN {threads} --genomeFastaFiles {input.genomeFasta} --genomeDir {params.outdir} --sjdbGTFfile {input.genomeAnno} 1>{log.log_file} 2>&1;
        '''

rule transcriptome_align_by_star:
    message:
        "Running STAR Alignment" 
    input:
        read1 = join(qc_outdir, qc_method, "{sample}", "{sample}.cleanR1.fq.gz"),
        read2 = join(qc_outdir, qc_method,"{sample}", "{sample}.cleanR2.fq.gz"),
        index_ok = join(junction_align_outdir, 'star', 'index', 'index.ok')
    output:
        sorted_bam = join(junction_align_outdir, 'star', 'align', '{sample}', '{sample}.sorted.bam'),
        counts = join(junction_align_outdir, 'star', 'align', '{sample}', '{sample}.ReadsPerGene.out.tab'),
    params:
        stranded=rna_library,
        index_outdir = join(junction_align_outdir, 'star', 'index'),
        mapping_outdir = join(junction_align_outdir, 'star', 'align', '{sample}'),
    threads: 8
    log:
        log_file = join(junction_align_outdir, 'star', 'align', '{sample}', "align.log.txt")
    priority: 10
    benchmark: join(junction_align_outdir, 'star', 'align', '{sample}', "benchmark.txt")
    shell:
        '''
        STAR --runMode alignReads --runThreadN {threads} --genomeDir {params.index_outdir} --readFilesIn {input.read1} {input.read2} --outFileNamePrefix {params.mapping_outdir}/{wildcards.sample}. --quantMode GeneCounts --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate 1>{log.log_file} 2>&1;
        ln -s {params.mapping_outdir}/{wildcards.sample}.Aligned.sortedByCoord.out.bam {output.sorted_bam} 1>>{log.log_file} 2>&1;
        '''

rule statistic_data_of_bam:
    message:
        '''
        ------------------------------
        after bam sorted, we need to statistic mapping rate or something
        ------------------------------
        '''
    input:
        sorted_bam = join(junction_align_outdir, junction_align_method, 'align', "{sample}", "{sample}.sorted.bam")
    output:
        bam_stat = join(junction_align_outdir, junction_align_method, 'stat', "{sample}", "{sample}.bam_stat.txt")
    benchmark:
        join(junction_align_outdir, junction_align_method, 'stat', "{sample}", "benchmark.txt")
    log:
        join(junction_align_outdir, junction_align_method, 'stat', "{sample}", "bam_stat.log")
    shell:
        '''
        bam_stat.py -i {input.sorted_bam} >{output.bam_stat} 2>{log};
        '''

rule mapping_report:
    message:
        '''
        ------------------------------
        mapping_report
        ------------------------------
        '''
    input:
        bam_stat = join(junction_align_outdir, junction_align_method, 'stat', "{sample}", "{sample}.bam_stat.txt")
    output:
        bam_report = report(join(config['reportsDir'], '4.mapping_stat', "{sample}.bam_stat.txt"), caption=join(snake_dir, "reports/cleanreads_stat.rst"), category="Step 3: clean reads mapping stat")
    shell:
        '''
        cp {input.bam_stat} {output.bam_report};
        '''

rule mapping:
    message:
        '''
        mapping result
        '''
    input:
        sorted_bam = expand(join(junction_align_outdir, junction_align_method, 'align', "{sample}", "{sample}.sorted.bam"), sample=SAMPLES),
        bam_report = expand(join(config['reportsDir'], '4.mapping_stat', "{sample}.bam_stat.txt"), sample=SAMPLES)
    output:
        mapping_ok = touch(join(flag_outdir, 'mapping.ok'))
