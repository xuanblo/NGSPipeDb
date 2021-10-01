sampling_method = config['sampling_method'] # tail, seqkit_number, seqkit_proportion, head, tail
sampling_data_outdir = join(config["resultsDir"], "sampling_reads")

rule sampling_data_by_links:
    message:
        '''
        ------------------------------
        use all of the reads, only create links to it need to be 
        ------------------------------
        '''
    input:
        read1 = join(config["samplesDir"], "{sample}"+config["read1Suffix"]),
        read2 = join(config["samplesDir"], "{sample}"+config["read2Suffix"])
    output:
        read1 = join(sampling_data_outdir, 'links', "{sample}", "{sample}"+config["read1Suffix"]),
        read2 = join(sampling_data_outdir, 'links', "{sample}", "{sample}"+config["read2Suffix"])
    threads:
        5
    benchmark:
        join(sampling_data_outdir, 'links', "{sample}", "benchmark.txt")
    log:
        join(sampling_data_outdir, 'links', "{sample}", "run.log.txt")
    run:
        import os
        read1_src = os.path.abspath(input.read1)
        read2_src = os.path.abspath(input.read2)
        os.symlink(read1_src, output.read1)
        os.symlink(read2_src, output.read2)

rule sampling_data_by_head:
    message:
        '''
        ------------------------------
        sampling raw reads by using first {} lines
        ------------------------------
        '''.format(config['sampling_value'])
    input:
        read1 = join(config["samplesDir"], "{sample}"+config["read1Suffix"]),
        read2 = join(config["samplesDir"], "{sample}"+config["read2Suffix"])
    output:
        read1 = join(sampling_data_outdir, 'head', "{sample}", "{sample}"+config["read1Suffix"]),
        read2 = join(sampling_data_outdir, 'head', "{sample}", "{sample}"+config["read2Suffix"])
    threads:
        5
    benchmark:
        join(sampling_data_outdir, 'head', "{sample}", "benchmark.txt")
    log:
        join(sampling_data_outdir, 'head', "{sample}", "run.log.txt")
    params:
        lines = config['sampling_value']
    shell:
        '''
        gunzip -c {input.read1} | head -{params.lines} | gzip > {output.read1} || sleep 0.1 2>{log};
        gunzip -c {input.read2} | head -{params.lines} | gzip > {output.read2} || sleep 0.1 2>>{log};
        '''

rule sampling_data_by_tail:
    message:
        '''
        ------------------------------
        sampling raw reads by using last {} lines
        ------------------------------
        '''.format(config['sampling_value'])
    input:
        read1 = join(config["samplesDir"], "{sample}"+config["read1Suffix"]),
        read2 = join(config["samplesDir"], "{sample}"+config["read2Suffix"])
    output:
        read1 = join(sampling_data_outdir, 'tail', "{sample}", "{sample}"+config["read1Suffix"]),
        read2 = join(sampling_data_outdir, 'tail', "{sample}", "{sample}"+config["read2Suffix"])
    threads:
        5
    benchmark:
        join(sampling_data_outdir, 'tail', "{sample}", "benchmark.txt")
    log:
        join(sampling_data_outdir, 'tail', "{sample}", "run.log.txt")
    params:
        lines = config['sampling_value']
    shell:
        '''
        gunzip -c {input.read1} | tail -{params.lines} | gzip > {output.read1} || sleep 0.1 2>{log};
        gunzip -c {input.read2} | tail -{params.lines} | gzip > {output.read2} || sleep 0.1 2>>{log};
        '''

rule sampling_data_by_seqkit_number:
    message:
        '''
        ------------------------------
        sampling raw reads by using {} reads, number of reads will be selected randomly 
        ------------------------------
        '''.format(config['sampling_value'])
    input:
        read1 = join(config["samplesDir"], "{sample}"+config["read1Suffix"]),
        read2 = join(config["samplesDir"], "{sample}"+config["read2Suffix"])
    output:
        read1 = join(sampling_data_outdir, 'seqkit_number', "{sample}", "{sample}"+config["read1Suffix"]),
        read2 = join(sampling_data_outdir, 'seqkit_number', "{sample}", "{sample}"+config["read2Suffix"])
    threads:
        5
    benchmark:
        join(sampling_data_outdir, 'seqkit_number', "{sample}", "benchmark.txt")
    log:
        join(sampling_data_outdir, 'seqkit_number', "{sample}", "run.log.txt")
    conda:
        "envs/seqkit.yaml"
    params:
        number = config['sampling_value']
    shell:
        '''
        gunzip -c {input.read1} | seqkit sample -n {params.number} -j {threads} -o {output.read1} 1>{log} 2>&1;
        gunzip -c {input.read2} | seqkit sample -n {params.number} -j {threads} -o {output.read2} 1>>{log} 2>&1;
        '''

rule sampling_data_by_seqkit_proportion:
    message:
        '''
        ------------------------------
        sampling raw reads by using {} proportion, reads will be selected randomly 
        ------------------------------
        '''.format(config['sampling_value'])
    input:
        read1 = join(config["samplesDir"], "{sample}"+config["read1Suffix"]),
        read2 = join(config["samplesDir"], "{sample}"+config["read2Suffix"])
    output:
        read1 = join(sampling_data_outdir, 'seqkit_proportion', "{sample}", "{sample}"+config["read1Suffix"]),
        read2 = join(sampling_data_outdir, 'seqkit_proportion', "{sample}", "{sample}"+config["read2Suffix"])
    threads:
        5
    benchmark:
        join(sampling_data_outdir, 'seqkit_proportion', "{sample}", "benchmark.txt")
    log:
        join(sampling_data_outdir, 'seqkit_proportion', "{sample}", "run.log.txt")
    conda:
        "envs/seqkit.yaml"
    params:
        proportion = config['sampling_value']
    shell:
        '''
        gunzip -c {input.read1} | seqkit sample -p {params.proportion} -j {threads} -o {output.read1} 1>{log} 2>&1;
        gunzip -c {input.read2} | seqkit sample -p {params.proportion} -j {threads} -o {output.read2} 1>>{log} 2>&1;
        '''

rule rawReads_data_stat:
    message:
        '''
        ------------------------------
        reads qc cleak mainly
        ------------------------------
        '''
    input:
        read1 = join(sampling_data_outdir, sampling_method, "{sample}", "{sample}"+config["read1Suffix"]),
        read2 = join(sampling_data_outdir, sampling_method, "{sample}", "{sample}"+config["read2Suffix"])
    output:
        read1 = join(sampling_data_outdir, 'stat', "{sample}", "{sample}_1.seqkit.summary"),
        read2 = join(sampling_data_outdir, 'stat', "{sample}", "{sample}_2.seqkit.summary"),
        read1_gc_qual = join(sampling_data_outdir, 'stat', "{sample}", "{sample}_1.seqkit.gc_qual"),
        read2_gc_qual = join(sampling_data_outdir, 'stat', "{sample}", "{sample}_2.seqkit.gc_qual")
    threads: 1
    benchmark:
        join(sampling_data_outdir, 'stat', "{sample}", "benchmark.txt")
    log:
        join(sampling_data_outdir, 'stat', "{sample}", "run.log.txt"),
    shell:
        '''
        seqkit stat -j {threads} -a {input.read1} >{output.read1} 2>{log};
        seqkit stat -j {threads} -a {input.read2} >{output.read2} 2>>{log};
        seqkit fx2tab -n -g -q -j {threads} {input.read1}|cut -f4-|awk '{{gc=gc+$1;qual=qual+$2;}}END{{print gc/NR"\\t"qual/NR}}' >{output.read1_gc_qual} 2>>{log};
        seqkit fx2tab -n -g -q -j {threads} {input.read2}|cut -f4-|awk '{{gc=gc+$1;qual=qual+$2;}}END{{print gc/NR"\\t"qual/NR}}' >{output.read2_gc_qual} 2>>{log};
        '''

rule rawReads_stat_merge:
    message:
        '''
        ------------------------------
        reads qc merge result
        ------------------------------
        '''
    input:
        read1 = expand(join(sampling_data_outdir, 'stat', "{sample}", "{sample}_1.seqkit.summary"), sample=SAMPLES),
        read2 = expand(join(sampling_data_outdir, 'stat', "{sample}", "{sample}_2.seqkit.summary"), sample=SAMPLES),
        read1_gc_qual = expand(join(sampling_data_outdir, 'stat', "{sample}", "{sample}_1.seqkit.gc_qual"), sample=SAMPLES),
        read2_gc_qual = expand(join(sampling_data_outdir, 'stat', "{sample}", "{sample}_2.seqkit.gc_qual"), sample=SAMPLES)
    output:
        rawdata_report = join(sampling_data_outdir, 'stat', "reads_product.csv"),
        rawReads_ok = touch(join(sampling_data_outdir, 'stat', 'statistic.completed'))
    benchmark:
        join(sampling_data_outdir, 'stat', "benchmark.txt")
    log:
        join(sampling_data_outdir, 'stat', "run.log.txt"),
    shell:
        '''
        python {snake_dir}/scripts/reads_product_stat.py {input.read1} {input.read2} {input.read1_gc_qual} {input.read2_gc_qual} {output.rawdata_report} 1>{log} 2>&1;
        '''

rule rawreads_stat_report:
    message:
        '''
        ------------------------------
        rawreads_stat_report
        ------------------------------
        '''
    input:
        rawdata_report = join(sampling_data_outdir, 'stat', "reads_product.csv")
    output:
        rawdata_report = report(join(config['reportsDir'], '2.rawreads_stat', "rawreads_product.csv"), caption=join(snake_dir, "reports/rawreads_stat.rst"), category="Step 1: rawreads stat & fastqc")
    shell:
        '''
        cp {input.rawdata_report} {output.rawdata_report};
        '''

rule sampling_reads:
    message:
        '''
        touch sampling_reads.ok
        '''
    input:
        read1 = expand(join(sampling_data_outdir, sampling_method, "{sample}", "{sample}"+config["read1Suffix"]), sample=SAMPLES),
        read2 = expand(join(sampling_data_outdir, sampling_method, "{sample}", "{sample}"+config["read2Suffix"]), sample=SAMPLES),
        stat = join(config['reportsDir'], '2.rawreads_stat', "rawreads_product.csv"),
    output:
        sampling_reads_ok = touch(join(flag_outdir, 'sampling_reads.ok'))
