# global description
report: join(snake_dir, "reports/workflow.rst")

"""
rule dag:
    input:
    output:
        "dag.svg"
    shell:
        '''
        snakemake --snakefile NGSPipeCode/Snakefile --configfile NGSPipeCode/config.yaml --dag | dot -Tsvg > dag.svg
        '''

rule pipeline:
    input:
        'dag.svg'
    output:
        report(join(config['reportsDir'], 'report', "fig1.svg"), caption="reports/fig1.rst", category="Step 1: 设计自己的pipeline")
    shell:
        '''
        cp {input} {output};
        '''
"""

#ruleorder: differential_expression_analysis_by_deseq2 > report_all

rule report_all:
    message:
        '''
        ------------------------------
        merge report result
        ------------------------------
        '''
    input:
        # 1. workflow
        workflow_report     = join(config['reportsDir'], '1.pipeline', "workflow.png"),
        # 2. rawdata statistic
        rawdata_report      = join(config['reportsDir'], '2.rawreads_stat', "rawreads_product.csv"),
        # 3. cleandata statistic
        cleandata_report    = join(config['reportsDir'], '3.cleanreads_stat', "cleanreads_product.csv"),
        # 4. mapping statistic
        bam_report          = expand(join(config['reportsDir'], '4.mapping_stat', "{sample}.bam_stat.txt"), sample=SAMPLES),
        # 5. expression statistic
        gene_pca_pdf        = join(config['reportsDir'], '5.exp_stat', "gene_expression_pca.pdf"),
        gene_heatmap_pdf    = join(config['reportsDir'], '5.exp_stat', "gene_sample_heatmap.pdf"),
        gene_sig_heatmap_pdf = join(config['reportsDir'], '5.exp_stat', "gene_sig_variance_heatmap.pdf"),
        #volcano_out_pdf = expand(join(config['reportsDir'], '5.exp_stat', "volcano", "{compare}"), compare=os.listdir(join(diff_outdir, "volcano"))),
        #volcano_out_pdf = join(config['reportsDir'], '5.exp_stat', "volcano"),
    output:
        touch(join(config['reportsDir'], "report.ok"))
    shell:
        '''

        '''

# 1. workflow
rule workflow_report:
    message:
        '''
        ------------------------------
        workflow_report
        ------------------------------
        '''
    input:
        workflow_report = join(snake_dir, "imgs", "workflow.png")
    output:
        pipeline_report = join(config['reportsDir'], '1.pipeline', "workflow.png")
    shell:
        '''
        if [ -e "{working_dir}/dag.png" ]; then
            mv {working_dir}/dag.png {output.pipeline_report};
        else
            cp {input.workflow_report} {output.pipeline_report};
        fi
        '''

# 2. rawdata statistic
rule rawreads_stat_report:
    message:
        '''
        ------------------------------
        rawreads_stat_report
        ------------------------------
        '''
    input:
        rawdata_report = join(rawReads_outdir, "reads_product.csv")
    output:
        rawdata_report = report(join(config['reportsDir'], '2.rawreads_stat', "rawreads_product.csv"), caption=join(snake_dir, "reports/rawreads_stat.rst"), category="Step 1: rawreads stat & fastqc")
    shell:
        '''
        cp {input.rawdata_report} {output.rawdata_report};
        '''

# 3. cleandata statistic
rule cleanreads_stat_report:
    message:
        '''
        ------------------------------
        cleanreads_stat_report (after qc)
        ------------------------------
        '''
    input:
        cleandata_report = join(cleanReads_outdir, "reads_product.csv")
    output:
        cleandata_report = report(join(config['reportsDir'], '3.cleanreads_stat', "cleanreads_product.csv"), caption=join(snake_dir, "reports/cleanreads_stat.rst"), category="Step 2: clean reads stat & fastqc")
    shell:
        '''
        cp {input.cleandata_report} {output.cleandata_report};
        '''

# 4. mapping statistic
rule cleanreads_mapping_report:
    message:
        '''
        ------------------------------
        cleanreads_mapping_report
        ------------------------------
        '''
    input:
        bam_stat = join(bam_outdir, "{sample}", "{sample}.bam_stat.txt")
    output:
        bam_report = report(join(config['reportsDir'], '4.mapping_stat', "{sample}.bam_stat.txt"), caption=join(snake_dir, "reports/cleanreads_stat.rst"), category="Step 3: clean reads mapping stat")
    shell:
        '''
        cp {input.bam_stat} {output.bam_report};
        '''

# 5. expression statistic

rule gene_expression_pca:
    message:
        '''
        ------------------------------
        gene_expression_pca
        ------------------------------
        '''
    input:
        pca_pdf = join(diff_outdir, "pca.pdf")
    output:
        gene_pca_pdf = report(join(config['reportsDir'], '5.exp_stat', "gene_expression_pca.pdf"), caption=join(snake_dir, "reports/expression.rst"), category="pca plot")
    shell:
        '''
        cp {input.pca_pdf} {output.gene_pca_pdf};
        '''

rule gene_expression_sample_heatmap:
    message:
        '''
        ------------------------------
        gene_expression_sample_heatmap
        ------------------------------
        '''
    input:
        heatmap_pdf = join(diff_outdir, "sample_heatmap.pdf")
    output:
        gene_heatmap_pdf = report(join(config['reportsDir'], '5.exp_stat', "gene_sample_heatmap.pdf"), caption=join(snake_dir, "reports/expression.rst"), category="sample heatmap")
    shell:
        '''
        cp {input.heatmap_pdf} {output.gene_heatmap_pdf};
        '''

rule top_gene_expression_variance_heatmap:
    message:
        '''
        ------------------------------
        top_gene_expression_variance_heatmap
        ------------------------------
        '''
    input:
        sig_heatmap_pdf = join(diff_outdir, "sig_gene_exp.pdf")
    output:
        gene_sig_heatmap_pdf = report(join(config['reportsDir'], '5.exp_stat', "gene_sig_variance_heatmap.pdf"), caption=join(snake_dir, "reports/expression.rst"), category="top_gene_expression_variance_heatmap")
    shell:
        '''
        cp {input.sig_heatmap_pdf} {output.gene_sig_heatmap_pdf};
        '''

rule differential_compare_volcano_plot:
    message:
        '''
        ------------------------------
        differential_compare_volcano_plot
        ------------------------------
        '''
    input:
        #volcano_pdf = join(diff_outdir, "volcano", "{compare}"),
        volcano_dir = join(diff_outdir, "volcano")
    output:
        #volcano_out_pdf = report(join(config['reportsDir'], '5.exp_stat', "volcano", "{compare}"), caption=join(snake_dir, "reports/cleanreads_stat.rst"), category="volcano")
        volcano_out_dir = report(join(config['reportsDir'], '5.exp_stat', "volcano"), caption=join(snake_dir, "reports/expression.rst"), category="volcano")
    shell:
        '''
        #cp {input.volcano_pdf} {output.volcano_out_pdf};
        cp -r {input.volcano_dir} {output.volcano_out_dir};
        '''

"""
rule assembled_transcript_report:
    message:
        '''
        ------------------------------
        assembled_transcript_report
        ------------------------------
        '''
    input:
        'NGSPipeCode/reports/fig2.png'
    output:
        report(join(config['reportsDir'], 'report', "fig2.png"), caption=join(snake_dir, "reports/fig2.rst"), category="Step 2: 统计")
    shell:
        "cp {input} {output}"


"""