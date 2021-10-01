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
        condition = config["conditionPath"]
    output:
        diff_outputok = touch(join(diff_outdir, "diff.ok")),
        pca_pdf = join(diff_outdir, 'pca.pdf'),
        sample_heatmap_pdf = join(diff_outdir, 'sample_heatmap.pdf'),
        sig_gene_exp_pdf = join(diff_outdir, 'sig_gene_exp.pdf'),
        volcano_dir = directory(join(diff_outdir, 'volcano')),
        normalized_counts = join(diff_outdir, 'normalized.counts.csv'),
    conda:
        '../envs/requirements_exp_r_env.yaml'
    benchmark:
        join(diff_outdir, "benchmark.txt")
    log:
        join(diff_outdir, "diff.log")
    shell:
        '''
        # run deseq2
        Rscript {snake_dir}/scripts/deseq2_call_diff_with_counts_matrix.R --countsMatrix {input.gene_count} --conditionFile {input.condition} --resultDir {diff_outdir} 1>{log} 2>&1;
        # draw pca
        Rscript {snake_dir}/scripts/plot_deseq2_sample_pca.R --deseq2RData {diff_outdir}/deseq2.RData --outpdf {diff_outdir}/pca.pdf 1>>{log} 2>&1;
        # draw sample heatmap
        Rscript {snake_dir}/scripts/plot_deseq2_sample_heatmap.R --deseq2RData {diff_outdir}/deseq2.RData --outpdf {diff_outdir}/sample_heatmap.pdf 1>>{log} 2>&1;
        # draw top variance gene exp heatmap
        Rscript {snake_dir}/scripts/plot_deseq2_sig_genes_heatmap.R --deseq2RData {diff_outdir}/deseq2.RData --outpdf {diff_outdir}/sig_gene_exp.pdf 1>>{log} 2>&1;
        # draw volcano plot
        Rscript {snake_dir}/scripts/plot_deseq2_enhance_volcano.R --deseq2RData {diff_outdir}/deseq2.RData --outdir {diff_outdir}/volcano 1>>{log} 2>&1;
        '''

if not "genomeAnno" in config.keys():
    config["genomeAnno"] = ''

rule diff_gene_add_note:
    message:
        '''
        ------------------------------
        add annotation to differential gene
        ------------------------------
        '''
    input:
        diff_outputok = join(diff_outdir, "diff.ok"),
        gtf = config["genomeAnno"]
    output:
        diff_outputok = touch(join(diff_anno_outdir, "diff.ok"))
    conda:
        '../envs/requirements_exp_r_env.yaml'
    benchmark:
        join(diff_anno_outdir, "benchmark.txt")
    log:
        join(diff_anno_outdir, "diff.log")
    shell:
        '''
        python {snake_dir}/scripts/add_anno2diff_result.py {input.gtf} {diff_outdir} {diff_anno_outdir} 1>{log} 2>&1;
        '''
