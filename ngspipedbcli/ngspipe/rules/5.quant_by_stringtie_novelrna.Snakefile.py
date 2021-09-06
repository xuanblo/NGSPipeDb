# ----------------------------------------------------------------------------
# restringtie
# ----------------------------------------------------------------------------
rule reStringtie_by_stringtie:
    message:
        '''
        stringtie was used to assemble transcript
        '''
    input:
        mergedGtf = join(transcript_assembly_outdir, "merged.gtf"),
        sorted_bam = join(junction_align_outdir, "{sample}", "{sample}.sorted.bam"),
    output:
        gtf = join(quantify_outdir, "{sample}", "{sample}.gtf"),
        expr = join(quantify_outdir, "{sample}", "{sample}.tab"),
        gene_count = join(quantify_outdir, "{sample}", "{sample}_gene_count_piece.csv")
    threads:
        4
    benchmark:
        join(quantify_outdir, "{sample}", "benchmark.txt")
    log:
        join(quantify_outdir, "{sample}", "run.log")
    shell:
        '''
        stringtie -p {threads} -G {input.mergedGtf} -o {output.gtf} -A {output.expr} -B -e -l {wildcards.sample} {input.sorted_bam} 1>{log} 2>&1;
        echo -e {wildcards.sample}\"\\t\"{output.gtf} >{output.gene_count} 2>>{log};
        '''

# ----------------------------------------------------------------------------
# counts
# ----------------------------------------------------------------------------
rule countMatrix_by_stringtie:
    message:
        '''
        fpkm, tpm, counts
        '''
    input:
        gene_count = expand(join(quantify_outdir, "{sample}", "{sample}_gene_count_piece.csv"), sample = SAMPLES),
        expr = expand(join(quantify_outdir, "{sample}", "{sample}.tab"), sample = SAMPLES)
    output:
        gene_count = join(quantify_outdir, "gene.csv"),
        transcript_count = join(quantify_outdir, "transcript.csv"),
        transcript_fpkm = join(quantify_outdir, "transcript_fpkm_all_samples.tsv"),
        gene_fpkm = join(quantify_outdir, "gene_fpkm_all_samples.tsv")
    log:
        join(quantify_outdir, "prepDE.log")
    benchmark:
        join(quantify_outdir, "benchmark.txt")
    shell:
        '''
        cat {input.gene_count} >{quantify_outdir}/gtflist.txt 2>{log};
        prepDE.py -i {quantify_outdir}/gtflist.txt -t {output.transcript_count} -g {output.gene_count} -l 150 >>{log} 2>&1;
        # TPM
        python {snake_dir}/scripts/stringtie_expression_matrix.py -m tpm -s {quantify_outdir} -t {quantify_outdir}/transcript_tpm_all_samples.tsv -g {quantify_outdir}/gene_tpm_all_samples.tsv >>{log} 2>&1;
        # FPKM
        python {snake_dir}/scripts/stringtie_expression_matrix.py -m fpkm -s {quantify_outdir} -t {quantify_outdir}/transcript_fpkm_all_samples.tsv -g {quantify_outdir}/gene_fpkm_all_samples.tsv >>{log} 2>&1;
        # Normalizition
        '''
