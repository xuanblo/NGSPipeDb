# ----------------------------------------------------------------------------
# refer_library_by_infer_experiment_RseQC
refer_library_outdir = join(config["result_dir"], "stat.rseqc", "refer_library_by_infer_experiment_RseQC")
hisat2_ourdir = join(config["result_dir"], "align.hisat2_stringtie", "step2_mapping_genome_by_hisat2")
# ----------------------------------------------------------------------------
rule refer_library_by_infer_experiment_RseQC:
    message:
        '''
        ------------------------------
        infer_experiment.py
        ------------------------------
        '''
    input:
        bam = join(hisat2_ourdir, "{sample}", "{sample}.sorted.bam"),
        genome_bed = config["genome_bed"]
    output:
        join(refer_library_outdir, "{sample}", "{sample}.strand.txt")
    params:
        dbPrefix = join(step1_outdir, step1_dbname),
        outdir = step1_outdir
    conda:
        py3env
    benchmark:
        join(refer_library_outdir, "{sample}", "benchmark.txt")
    log:
        join(refer_library_outdir, "{sample}", "infer.log")
    shell:
        '''
        infer_experiment.py -r {input.genome_bed} -i {input.bam} > {output}
        '''

# ----------------------------------------------------------------------------
# step1_bam_stat_outdir
step1_bam_stat_outdir = join(config["result_dir"], "stat.rseqc", "bam_stat_by_bam_stat_RseQC")
hisat2_ourdir = join(config["result_dir"], "align.hisat2_stringtie", "step2_mapping_genome_by_hisat2")
# ----------------------------------------------------------------------------
rule step1_bam_stat_by_bam_stat_RseQC:
    message:
        '''
        ------------------------------
        bam_stat.py
        ------------------------------
        '''
    input:
        join(hisat2_ourdir, "{sample}", "{sample}.sorted.bam"),
    output:
        join(step1_bam_stat_outdir, "{sample}", "{sample}.bam_stat.txt")
    conda:
        py3env
    benchmark:
        join(step1_bam_stat_outdir, "{sample}", "benchmark.txt")
    log:
        join(step1_bam_stat_outdir, "{sample}", "bam_stat.log")
    shell:
        '''
        bam_stat.py -i {input} >{output} 2>{log};
        '''


# ----------------------------------------------------------------------------
# step2_bam_stat_outdir
step2_bam_stat_outdir = join(config["result_dir"], "stat.rseqc", "bam_stat_by_bam_stat_RseQC")
# ----------------------------------------------------------------------------
rule step2_bam_stat_merge_by_script:
    message:
        '''
        ------------------------------
        script/bam_stat_merge.py
        ------------------------------
        '''
    input:
        expand(join(step1_bam_stat_outdir, "{sample}", "{sample}.bam_stat.txt"), sample=SAMPLES)
    output:
        join(step2_bam_stat_outdir, "bam_stat_merge.csv")
    conda:
        py3env
    benchmark:
        join(step2_bam_stat_outdir, "benchmark.txt")
    log:
        join(step2_bam_stat_outdir, "bam_stat_merge.log")
    shell:
        '''
        python script/bam_stat_merge.py {input} {output} 2>{log};
        '''
