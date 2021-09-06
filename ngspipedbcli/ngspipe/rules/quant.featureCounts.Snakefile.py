# ----------------------------------------------------------------------------
# step1
featureCounts_outdir = join(config["result_dir"], "quant.featureCounts", "feature_quant_by_featureCounts")
hisat2_ourdir = join(config["result_dir"], "align.hisat2_stringtie", "step2_mapping_genome_by_hisat2")
if config["rna_library"] == "fr-firststrand":
    featureCounts_library = "-s 2"
elif config["rna_library"] == "fr-secondstrand":
    featureCounts_library = "-s 1"
elif config["rna_library"] == "none":
    featureCounts_library = ""
# ----------------------------------------------------------------------------
rule quant_counts_by_feature_counts:
    input:
        genome_gtf = config["genome_gtf"],
        bam = expand(join(hisat2_ourdir, "{sample}", "{sample}.sorted.bam"), sample=SAMPLES)
    output:
        gene_counts = join(featureCounts_outdir, "gene.counts.tsv"),
        transcript_counts = join(featureCounts_outdir, "transcript.counts.tsv"),
    params:
        featureCounts_library = featureCounts_library
    conda:
        py3env
    threads:
        5
    benchmark:
        join(featureCounts_outdir, "{sample}", "benchmark.txt")
    log:
        join(featureCounts_outdir, "{sample}", "feature_counts.log")
    shell:
        '''
        featureCounts -g gene_id -t exon -p -a {input.genome_gtf} -o {output.gene_counts} -s 2 -T {threads} {input.bam} 1>{log} 2>&1;
        featureCounts -g transcript_id -t exon -p -a {input.genome_gtf} -o {output.transcript_counts} -s 2 -T {threads} {input.bam} 1>>{log} 2>&1;
        '''

# ----------------------------------------------------------------------------
# step2
merge_counts_ourdir = join(config["result_dir"], "quant.featureCounts", "merge_featurecount_by_script")
# ----------------------------------------------------------------------------
rule merge_featurecount_by_script:
    input:
        expand(join(featureCounts_outdir, "{sample}", "{sample}.counts_stat.txt"), sample=SAMPLES)
        #rules.quant_counts_by_feature_counts.output
    output:
        join(merge_counts_ourdir, "featurecount.tsv")
    conda:
        py3env
    benchmark:
        join(merge_counts_ourdir, "benchmark.txt")
    log:
        join(merge_counts_ourdir, "feature_counts.log")
    shell:
        '''
        python script/merge_featureCount.py {input} {output} 1>{log} 2>&1;
        '''