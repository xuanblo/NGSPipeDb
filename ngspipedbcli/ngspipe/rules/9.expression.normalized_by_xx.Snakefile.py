# ----------------------------------------------------------------------------
# input are raw count, replicates must >=3
# ----------------------------------------------------------------------------
rule differential_expression_analysis_by_deseq2:
    message:
        '''
        ------------------------------
        differential_expression_analysis_by_deseq2
        ------------------------------
        '''
    input:
        gene_count = join(quantify_outdir, "gene.csv"),
        condition = config["condition"]
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
        Rscript ../scripts/differential_expression_use_deseq2.R --count_matrix {input.gene_count} --condition {input.condition} --outdir {} 1>{log} 2>&1;
        '''