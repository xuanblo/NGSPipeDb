rule samtools_remove_duplication:
    input:
        bwa_bam_sort = join(mapping_outdir, "{sample}", "{sample}.sorted.bam"),
    output:
        rmdup_bam = join(remove_duplication_outdir, "{sample}", "{sample}.sorted_rmdup.bam"),
    log:
        join(remove_duplication_outdir, "{sample}", "{sample}.log")
    threads: 4
    shell:
        '''
        samtools sort -n {input.bwa_bam_sort}|samtools fixmate -m - -|samtools sort -|samtools markdup -r - {output.rmdup_bam} 1>{log} 2>&1;
        '''