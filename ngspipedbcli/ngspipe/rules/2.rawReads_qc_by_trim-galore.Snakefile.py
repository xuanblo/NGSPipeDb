# ----------------------------------------------------------------------------
# rawReads_qc_by_trim_galore

# ----------------------------------------------------------------------------
rule rawReads_qc_by_trim_galore:
    message:
        '''
        ------------------------------
        we need to use clean reads for next step mapping
        ------------------------------
        '''
    input:
        read1 = join(sampling_data_outdir, "{sample}"+config["read1Suffix"]),
        read2 = join(sampling_data_outdir, "{sample}"+config["read2Suffix"])
    output:
        read1 = join(qc_outdir, "{sample}", "{sample}.cleanR1.fq.gz"),
        read2 = join(qc_outdir, "{sample}", "{sample}.cleanR2.fq.gz")
    log:
        join(qc_outdir, "{sample}", "qc.log")
    benchmark:
        join(qc_outdir, "{sample}", "benchmark.txt")
    threads: 5
    shell:
        '''
        trim_galore -j {threads} -q 30 --phred33 --gzip --length 36 --basename {wildcards.sample} --fastqc --paired {input.read1} {input.read2} -o {qc_outdir}/{wildcards.sample} 1>{log} 2>&1;
        ln -s {wildcards.sample}_val_1.fq.gz {output.read1} 1>>{log} 2>&1;
        ln -s {wildcards.sample}_val_2.fq.gz {output.read2} 1>>{log} 2>&1;
        '''


