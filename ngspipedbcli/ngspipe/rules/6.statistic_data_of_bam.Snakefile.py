# ----------------------------------------------------------------------------
# statistic_data_of_bam
# ----------------------------------------------------------------------------
rule statistic_data_of_bam:
    message:
        '''
        ------------------------------
        after bam sorted, we need to statistic mapping rate or something
        ------------------------------
        '''
    input:
        sorted_bam = join(junction_align_outdir, "{sample}", "{sample}.sorted.bam")
    output:
        bam_stat = join(bam_outdir, "{sample}", "{sample}.bam_stat.txt")
    benchmark:
        join(bam_outdir, "{sample}", "benchmark.txt")
    log:
        join(bam_outdir, "{sample}", "bam_stat.log")
    shell:
        '''
        bam_stat.py -i {input.sorted_bam} >{output.bam_stat} 2>{log};
        '''

rule statistic_data_of_bam_ok:
    input:
        bams_stat = expand(join(bam_outdir, "{sample}", "{sample}.bam_stat.txt"), sample=SAMPLES)
    output:
        bam_ok = touch(join(bam_outdir, 'statistic.completed'))
    shell:
        ""
