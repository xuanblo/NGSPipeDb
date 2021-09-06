# ----------------------------------------------------------------------------
# step1
htseq_quant_outdir = join(config["result_dir"], "quant.HtSeqCount", "feature_quant_by_htseq_counts")
hisat2_ourdir = join(config["result_dir"], "align.hisat2_stringtie", "step2_mapping_genome_by_hisat2")
if config["rna_library"] == "fr-firststrand":
    htseq_library = "-s reverse"
elif config["rna_library"] == "fr-secondstrand":
    htseq_library = "-s forward"
elif config["rna_library"] == "none":
    htseq_library = ""
# ----------------------------------------------------------------------------
rule feature_quant_step1_by_htseq_counts:
    input:
        genome_gtf = config["genome_gtf"],
        bam = join(hisat2_ourdir, "{sample}", "{sample}.sorted.bam")
    output:
        join(htseq_quant_outdir, "{sample}", "{sample}.htseqcount.txt")
    params:
        htseq_library = htseq_library
    conda:
        py3env
    threads:
        5
    benchmark:
        join(htseq_quant_outdir, "{sample}", "benchmark.txt")
    log:
        join(htseq_quant_outdir, "{sample}", "feature_quant.log")
    shell:
        '''
        htseq-count {params.htseq_library} -f bam {input.bam} {input.genome_gtf} >{output} 2>{log};
        '''

# ----------------------------------------------------------------------------
# step2
htseq_merge_counts_ourdir = join(config["result_dir"], "quant.HtSeqCount", "merge_counts_by_script")
# ----------------------------------------------------------------------------
rule merge_counts_by_script:
    input:
        expand(join(htseq_quant_outdir, "{sample}", "{sample}.htseqcount.txt"), sample=SAMPLES)
    output:
        join(htseq_merge_counts_ourdir, "htseqcount.tsv")
    conda:
        py3env
    benchmark:
        join(htseq_merge_counts_ourdir, "benchmark.txt")
    log:
        join(htseq_merge_counts_ourdir, "merge_counts.log")
    shell:
        '''
        python script/merge_htseqcount.py {input} {output} 1>{log} 2>&1;
        '''
    
rule htseqcount_report:
    input:
        counts_stat = expand(config["resultfolder"] + "/htseqcount/{sample}/{sample}.htseqcount.txt", sample=SAMPLES)
    output:
        htseqcount_report = config["resultfolder"] + "/htseqcount_stat_report/reads_counts_stat.csv"
    shell:
        "python script/reads_distrubution.py {input.counts_stat} {output.htseqcount_report}"
        
rule htseqcount:
    input:
        gff = config["mRNA"],
        bam = config["resultfolder"] + "/mapping/{sample}/{sample}.sorted.bam"
    output:
        htseqcount = config["resultfolder"] + "/htseqcount/{sample}/{sample}.htseqcount.txt"
    threads: 2
    log:
        config["resultfolder"] + "/htseqcount/{sample}/{sample}.log"
    shell:
        "htseq-count -s reverse -f bam {input.bam} {input.gff} >{output.htseqcount} 2>{log}"

