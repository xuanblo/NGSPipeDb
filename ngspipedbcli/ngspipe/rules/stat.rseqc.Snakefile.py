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

rule bam_stat:
    input:
        bam = config["resultfolder"] + "/mapping/{sample}/{sample}.sorted.bam",
        gff = config["gtf"],
        genome = config["genome"],
    output:
        bam_stat = config["resultfolder"] + "/mapping_stat/{sample}/{sample}.bam_stat.txt"
    log:
        config["resultfolder"] + "/mapping_stat/{sample}/{sample}.log"
    params:
        outdir = config["resultfolder"] + "/mapping_stat/{sample}"
    conda: "env/RSeQC.yaml"
    threads: 5
    shell:
        '''
        # RSeQC
        bam_stat.py -i {input.bam} >{output.bam_stat} 2>{log};
        # qualimap
        qualimap bamqc -bam {input.bam} --java-mem-size=4G --sequencing-protocol strand-specific-reverse -nt {threads} --collect-overlap-pairs -gff {input.gff} -outdir {params.outdir}/qualimap --output-genome-coverage {params.outdir}/qualimap/genome.coverage.txt
        '''

rule bam_stat_report:
    input:
        bams_stat = expand(config["resultfolder"] + "/mapping_stat/{sample}/{sample}.bam_stat.txt", sample=SAMPLES)
    output:
        bam_stat_report = config["resultfolder"] + "/stat_report/bam_stat_report/reads_mapping_stat.csv"
    log: config["resultfolder"] + "/stat_report/bam_stat_report/reads_mapping_stat.log"
    shell:
        "python script/reads_mapping_stat.py {input.bams_stat} {output.bam_stat_report}"