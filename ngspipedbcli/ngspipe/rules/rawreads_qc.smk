
qc_method = config['qc_method'] # trimmomatic or cutadate or fastp or fastqc or multiqc
qc_outdir = join(config["resultsDir"], "qc")

rule rawReads_qc_by_trim_galore:
    message:
        '''
        ------------------------------
        we need to use clean reads for next step mapping
        ------------------------------
        '''
    input:
        read1 = join(sampling_data_outdir, sampling_method, "{sample}", "{sample}"+config["read1Suffix"]),
        read2 = join(sampling_data_outdir, sampling_method, "{sample}", "{sample}"+config["read2Suffix"])
    output:
        read1 = join(qc_outdir, 'trim-galore', "{sample}", "{sample}.cleanR1.fq.gz"),
        read2 = join(qc_outdir, 'trim-galore',"{sample}", "{sample}.cleanR2.fq.gz")
    log:
        join(qc_outdir, 'trim-galore', "{sample}", "qc.log")
    benchmark:
        join(qc_outdir, 'trim-galore', "{sample}", "benchmark.txt")
    threads: 5
    shell:
        '''
        trim_galore -j {threads} -q 30 --phred33 --gzip --length 36 --basename {wildcards.sample} --fastqc --paired {input.read1} {input.read2} -o {qc_outdir}/trim-galore/{wildcards.sample} 1>{log} 2>&1;
        ln -s {wildcards.sample}_val_1.fq.gz {output.read1} 1>>{log} 2>&1;
        ln -s {wildcards.sample}_val_2.fq.gz {output.read2} 1>>{log} 2>&1;
        '''

rule rawReads_qc_by_trimmomatic:
    message:
        '''
        ------------------------------
        we need to use clean reads for next step mapping use trimmomatic
        ------------------------------
        '''
    input:
        read1 = join(sampling_data_outdir, sampling_method, "{sample}", "{sample}"+config["read1Suffix"]),
        read2 = join(sampling_data_outdir, sampling_method, "{sample}", "{sample}"+config["read2Suffix"])
    output:
        read1_paired = join(qc_outdir, 'trimmomatic', "{sample}", "{sample}.cleanR1.fq.gz"),
        read2_paired = join(qc_outdir, 'trimmomatic', "{sample}", "{sample}.cleanR2.fq.gz"),
        read1_unpaired = join(qc_outdir, 'trimmomatic', "{sample}", "{sample}.cleanR1.unpaired.fq.gz"),
        read2_unpaired = join(qc_outdir, 'trimmomatic', "{sample}", "{sample}.cleanR2.unpaired.fq.gz"),
    log:
        join(qc_outdir, 'trimmomatic', "{sample}", "qc.log.txt")
    benchmark:
        join(qc_outdir, 'trimmomatic', "{sample}", "benchmark.txt")
    threads: 5
    params:
        trimmer=["TRAILING:3"],
    shell:
        '''
        trimmomatic PE -threads {threads} {input.read1} {input.read2} {output.read1_paired} {output.read1_unpaired} {output.read2_paired} {output.read2_unpaired} {params.trimmer} 1>{log} 2>&1;
        '''

rule cleanReads_data_stat:
    message:
        '''
        ------------------------------
        reads qc cleak mainly
        ------------------------------
        '''
    input:
        read1 = join(qc_outdir, qc_method, "{sample}",  "{sample}.cleanR1.fq.gz"),
        read2 = join(qc_outdir, qc_method, "{sample}", "{sample}.cleanR2.fq.gz")
    output:
        read1 = join(qc_outdir, 'stat', "{sample}", "{sample}_1.seqkit.summary"),
        read2 = join(qc_outdir, 'stat', "{sample}", "{sample}_2.seqkit.summary"),
        read1_gc_qual = join(qc_outdir, 'stat', "{sample}", "{sample}_1.seqkit.gc_qual"),
        read2_gc_qual = join(qc_outdir, 'stat', "{sample}", "{sample}_2.seqkit.gc_qual")
    threads: 1
    benchmark:
        join(qc_outdir, 'stat', "{sample}", "benchmark.txt")
    log:
        join(qc_outdir, 'stat', "{sample}", "run.log.txt"),
    shell:
        '''
        seqkit stat -j {threads} -a {input.read1} >{output.read1} 2>{log};
        seqkit stat -j {threads} -a {input.read2} >{output.read2} 2>>{log};
        seqkit fx2tab -n -g -q -j {threads} {input.read1}|cut -f4-|awk '{{gc=gc+$1;qual=qual+$2;}}END{{print gc/NR"\\t"qual/NR}}' >{output.read1_gc_qual} 2>>{log};
        seqkit fx2tab -n -g -q -j {threads} {input.read2}|cut -f4-|awk '{{gc=gc+$1;qual=qual+$2;}}END{{print gc/NR"\\t"qual/NR}}' >{output.read2_gc_qual} 2>>{log};
        '''

rule cleanReads_stat_merge:
    message:
        '''
        ------------------------------
        reads qc merge result
        ------------------------------
        '''
    input:
        read1 = expand(join(qc_outdir, 'stat', "{sample}", "{sample}_1.seqkit.summary"), sample=SAMPLES),
        read2 = expand(join(qc_outdir, 'stat', "{sample}", "{sample}_2.seqkit.summary"), sample=SAMPLES),
        read1_gc_qual = expand(join(qc_outdir, 'stat', "{sample}", "{sample}_1.seqkit.gc_qual"), sample=SAMPLES),
        read2_gc_qual = expand(join(qc_outdir, 'stat', "{sample}", "{sample}_2.seqkit.gc_qual"), sample=SAMPLES)
    output:
        rawdata_report = join(qc_outdir, 'stat', "reads_product.csv"),
        rawReads_ok = touch(join(qc_outdir, 'stat', 'statistic.completed'))
    benchmark:
        join(qc_outdir, 'stat', "benchmark.txt")
    log:
        join(qc_outdir, 'stat', "run.log.txt"),
    shell:
        '''
        python {snake_dir}/scripts/reads_product_stat.py {input.read1} {input.read2} {input.read1_gc_qual} {input.read2_gc_qual} {output.rawdata_report} 1>{log} 2>&1;
        '''

rule cleanreads_stat_report:
    message:
        '''
        ------------------------------
        cleanreads_stat_report (after qc)
        ------------------------------
        '''
    input:
        cleandata_report = join(qc_outdir, 'stat', "reads_product.csv")
    output:
        cleandata_report = report(join(config['reportsDir'], '3.cleanreads_stat', "cleanreads_product.csv"), caption=join(snake_dir, "reports/cleanreads_stat.rst"), category="Step 2: cleanreads stat & fastqc")
    shell:
        '''
        cp {input.cleandata_report} {output.cleandata_report};
        '''

rule rawreads_qc:
    message:
        '''
        -------------------
        final clean reads
        -------------------
        '''
    input:
        read1 = expand(join(qc_outdir, qc_method, "{sample}", "{sample}.cleanR1.fq.gz"), sample=SAMPLES),
        read2 = expand(join(qc_outdir, qc_method, "{sample}", "{sample}.cleanR2.fq.gz"), sample=SAMPLES),
        cleandata_report = join(config['reportsDir'], '3.cleanreads_stat', "cleanreads_product.csv")
    output:
        rawreads_qc_ok = touch(join(flag_outdir, 'rawreads_qc.ok'))