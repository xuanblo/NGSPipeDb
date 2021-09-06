# ----------------------------------------------------------------------------
# step1 bowtie2 build

# ----------------------------------------------------------------------------
rule call_variation_by_bcftools:
    message:
        '''
        ------------------------------
        samtools-bcftools
        ------------------------------
        '''
    input:
        genomeFasta = config["genomeFasta"],
        sorted_bam = join(mapping_outdir, "{sample}", "{sample}.sorted.bam"),
    output:
        vcf = join(variation_outdir, "{sample}", "{sample}.vcf"),
    benchmark:
        join(variation_outdir, "{sample}", "benchmark.txt")
    log:
        join(variation_outdir, "{sample}", "call.log")
    threads:
        5
    shell:
        '''
        bcftools mpileup -Ou -f {input.genomeFasta} {input.sorted_bam} |bcftools call -Ou -mv | bcftools filter -s LowQual -e '%QUAL<20 || DP>100' >{output.vcf} 2>{log};
        '''

