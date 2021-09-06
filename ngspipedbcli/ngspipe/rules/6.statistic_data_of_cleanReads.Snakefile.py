rule cleanReads_data_stat:
    message:
        '''
        ------------------------------
        reads qc cleak mainly
        ------------------------------
        '''
    input:
        read1 = join(qc_outdir, "{sample}", "{sample}.cleanR1.fq.gz"),
        read2 = join(qc_outdir, "{sample}", "{sample}.cleanR2.fq.gz")
    output:
        read1 = join(cleanReads_outdir, "{sample}", "{sample}_1.seqkit.summary"),
        read2 = join(cleanReads_outdir, "{sample}", "{sample}_2.seqkit.summary"),
        read1_gc_qual = join(cleanReads_outdir, "{sample}", "{sample}_1.seqkit.gc_qual"),
        read2_gc_qual = join(cleanReads_outdir, "{sample}", "{sample}_2.seqkit.gc_qual")
    threads: 1
    benchmark:
        join(cleanReads_outdir, "{sample}", "benchmark.txt")
    log:
        join(cleanReads_outdir, "{sample}", "run.log"),
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
        read1 = expand(join(cleanReads_outdir, "{sample}", "{sample}_1.seqkit.summary"), sample=SAMPLES),
        read2 = expand(join(cleanReads_outdir, "{sample}", "{sample}_2.seqkit.summary"), sample=SAMPLES),
        read1_gc_qual = expand(join(cleanReads_outdir, "{sample}", "{sample}_1.seqkit.gc_qual"), sample=SAMPLES),
        read2_gc_qual = expand(join(cleanReads_outdir, "{sample}", "{sample}_2.seqkit.gc_qual"), sample=SAMPLES)
    output:
        cleandata_report = join(cleanReads_outdir, "reads_product.csv"),
        cleanReads_ok = touch(join(cleanReads_outdir, 'statistic.completed'))
    benchmark:
        join(cleanReads_outdir, "benchmark.txt")
    log:
        join(cleanReads_outdir, "run.log"),
    shell:
        '''
        python {snake_dir}/scripts/reads_product_stat.py {input.read1} {input.read2} {input.read1_gc_qual} {input.read2_gc_qual} {output.cleandata_report} 1>{log} 2>&1;
        '''