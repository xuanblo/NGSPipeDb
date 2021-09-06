rule rawReads_data_stat:
    message:
        '''
        ------------------------------
        reads qc cleak mainly
        ------------------------------
        '''
    input:
        read1 = join(sampling_data_outdir, "{sample}"+config["read1Suffix"]),
        read2 = join(sampling_data_outdir, "{sample}"+config["read2Suffix"])
    output:
        read1 = join(rawReads_outdir, "{sample}", "{sample}_1.seqkit.summary"),
        read2 = join(rawReads_outdir, "{sample}", "{sample}_2.seqkit.summary"),
        read1_gc_qual = join(rawReads_outdir, "{sample}", "{sample}_1.seqkit.gc_qual"),
        read2_gc_qual = join(rawReads_outdir, "{sample}", "{sample}_2.seqkit.gc_qual")
    threads: 1
    benchmark:
        join(rawReads_outdir, "{sample}", "benchmark.txt")
    log:
        join(rawReads_outdir, "{sample}", "run.log"),
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
        read1 = expand(join(rawReads_outdir, "{sample}", "{sample}_1.seqkit.summary"), sample=SAMPLES),
        read2 = expand(join(rawReads_outdir, "{sample}", "{sample}_2.seqkit.summary"), sample=SAMPLES),
        read1_gc_qual = expand(join(rawReads_outdir, "{sample}", "{sample}_1.seqkit.gc_qual"), sample=SAMPLES),
        read2_gc_qual = expand(join(rawReads_outdir, "{sample}", "{sample}_2.seqkit.gc_qual"), sample=SAMPLES)
    output:
        rawdata_report = join(rawReads_outdir, "reads_product.csv"),
        rawReads_ok = touch(join(rawReads_outdir, 'statistic.completed'))
    benchmark:
        join(rawReads_outdir, "benchmark.txt")
    log:
        join(rawReads_outdir, "run.log"),
    shell:
        '''
        python {snake_dir}/scripts/reads_product_stat.py {input.read1} {input.read2} {input.read1_gc_qual} {input.read2_gc_qual} {output.rawdata_report} 1>{log} 2>&1;
        '''