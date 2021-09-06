rule picard_remove_duplication:
    input:
        bwa_bam_sort = join(mapping_outdir, "{sample}", "{sample}.sorted.bam"),
    output:
        rmdup_bam = join(remove_duplication_outdir, "{sample}", "{sample}.sorted_rmdup.bam"),
        rmdup_matrix = join(remove_duplication_outdir, "{sample}", "{sample}.sorted_rmdup.matrix")
    log:
        join(remove_duplication_outdir, "{sample}", "{sample}.log")
    threads: 4
    shell:
        '''
        picard MarkDuplicates REMOVE_DUPLICATES=true USE_JDK_DEFLATER=true USE_JDK_INFLATER=true I={input.bwa_bam_sort} O={output.rmdup_bam} M={output.rmdup_matrix} 1>{log} 2>&1;
        '''