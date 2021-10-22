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

rule gene_expression_pca_report:
    message:
        '''
        ------------------------------
        gene_expression_pca
        ------------------------------
        '''
    input:
        pca_pdf = join(diff_outdir, 'DESeq2', "pca.pdf")
    output:
        gene_pca_pdf = report(join(config['reportsDir'], '5.exp_stat', "gene_expression_pca.pdf"), caption=join(snake_dir, "reports/expression.rst"), category="pca plot")
    log: join(config['reportsDir'], '5.exp_stat', "gene_expression_pca.log.txt")
    shell:
        '''
        cp {input.pca_pdf} {output.gene_pca_pdf} 1>{log} 2>&1;
        '''

rule gene_expression_sample_heatmap_report:
    message:
        '''
        ------------------------------
        gene_expression_sample_heatmap
        ------------------------------
        '''
    input:
        heatmap_pdf = join(diff_outdir, 'DESeq2', "sample_heatmap.pdf")
    output:
        gene_heatmap_pdf = report(join(config['reportsDir'], '5.exp_stat', "gene_sample_heatmap.pdf"), caption=join(snake_dir, "reports/expression.rst"), category="sample heatmap")
    log: join(config['reportsDir'], '5.exp_stat', "gene_sample_heatmap.log.txt")
    shell:
        '''
        cp {input.heatmap_pdf} {output.gene_heatmap_pdf} 1>{log} 2>&1;
        '''

rule top_gene_expression_variance_heatmap_report:
    message:
        '''
        ------------------------------
        top_gene_expression_variance_heatmap
        ------------------------------
        '''
    input:
        sig_heatmap_pdf = join(diff_outdir, 'DESeq2', "sig_gene_exp.pdf")
    output:
        gene_sig_heatmap_pdf = report(join(config['reportsDir'], '5.exp_stat', "gene_sig_variance_heatmap.pdf"), caption=join(snake_dir, "reports/expression.rst"), category="top_gene_expression_variance_heatmap")
    log: join(config['reportsDir'], '5.exp_stat', "gene_sig_variance_heatmap.log.txt")
    shell:
        '''
        cp {input.sig_heatmap_pdf} {output.gene_sig_heatmap_pdf} 1>{log} 2>&1;
        '''

rule differential_compare_volcano_plot_report:
    message:
        '''
        ------------------------------
        differential_compare_volcano_plot
        ------------------------------
        '''
    input:
        volcano_dir = join(diff_outdir, 'DESeq2', "volcano"),
    output:
        report(directory(join(config['reportsDir'], '5.exp_stat', "volcano")), patterns=["{name}.volcano.pdf"], caption=join(snake_dir, "reports/expression.rst"), category="volcano"),
    shell:
        '''
        mkdir -p {output};
        for i in `ls {input.volcano_dir}|sed 's/.volcano.pdf//'`;do cp {input.volcano_dir}/$i.volcano.pdf {output}/$i.volcano.pdf;done
        '''

rule differential_expression:
    message:
        '''
        differential_expression ok
        '''
    input:
        quantification_ok = join(flag_outdir, 'quantification.ok'),
        normalized_counts = join(diff_outdir, diff_method, 'normalized.counts.csv'),
        gene_pca_pdf = join(config['reportsDir'], '5.exp_stat', "gene_expression_pca.pdf"),
        gene_heatmap_pdf = join(config['reportsDir'], '5.exp_stat', "gene_sample_heatmap.pdf"),
        gene_sig_heatmap_pdf = join(config['reportsDir'], '5.exp_stat', "gene_sig_variance_heatmap.pdf"),
        volcano_out_pdf = join(config['reportsDir'], '5.exp_stat', "volcano"),
    output:
        differential_expression_ok = touch(join(flag_outdir, 'differential_expression.ok'))
