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
