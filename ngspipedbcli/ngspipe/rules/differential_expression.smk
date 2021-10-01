diff_method = config['diff_method']
diff_outdir = join(config["resultsDir"], "diff")

rule differential_expression_analysis_by_deseq2:
    message:
        '''
        ------------------------------
        differential_expression_analysis_by_deseq2
        ------------------------------
        '''
    input:
        gene_count = join(quantify_outdir, quantify_method, 'quant', "gene.counts.csv"),
        condition = new_conditionpath
    output:
        diff_outputok = touch(join(diff_outdir, 'DESeq2', "diff.ok")),
        pca_pdf = join(diff_outdir, 'DESeq2', 'pca.pdf'),
        sample_heatmap_pdf = join(diff_outdir, 'DESeq2', 'sample_heatmap.pdf'),
        sig_gene_exp_pdf = join(diff_outdir, 'DESeq2', 'sig_gene_exp.pdf'),
        volcano_dir = directory(join(diff_outdir, 'DESeq2', 'volcano')),
        normalized_counts = join(diff_outdir, 'DESeq2', 'normalized.counts.csv'),
    benchmark:
        join(diff_outdir, 'DESeq2', "benchmark.txt")
    log:
        join(diff_outdir, 'DESeq2', "diff.log.txt")
    shell:
        '''
        # run deseq2
        Rscript {snake_dir}/scripts/deseq2_call_diff_with_counts_matrix.R --countsMatrix {input.gene_count} --conditionFile {input.condition} --resultDir {diff_outdir}/DESeq2 1>{log} 2>&1;
        # draw pca
        Rscript {snake_dir}/scripts/plot_deseq2_sample_pca.R --deseq2RData {diff_outdir}/DESeq2/deseq2.RData --outpdf {output.pca_pdf} 1>>{log} 2>&1;
        # draw sample heatmap
        Rscript {snake_dir}/scripts/plot_deseq2_sample_heatmap.R --deseq2RData {diff_outdir}/DESeq2/deseq2.RData --outpdf {output.sample_heatmap_pdf} 1>>{log} 2>&1;
        # draw top variance gene exp heatmap
        Rscript {snake_dir}/scripts/plot_deseq2_sig_genes_heatmap.R --deseq2RData {diff_outdir}/DESeq2/deseq2.RData --outpdf {output.sig_gene_exp_pdf} 1>>{log} 2>&1;
        # draw volcano plot
        Rscript {snake_dir}/scripts/plot_deseq2_enhance_volcano.R --deseq2RData {diff_outdir}/DESeq2/deseq2.RData --outdir {output.volcano_dir} 1>>{log} 2>&1;
        '''



rule differential_expression:
    message:
        '''
        differential_expression ok
        '''
    input:
        normalized_counts = join(diff_outdir, diff_method, 'normalized.counts.csv'),
    output:
        differential_expression_ok = touch(join(flag_outdir, 'differential_expression.ok'))
