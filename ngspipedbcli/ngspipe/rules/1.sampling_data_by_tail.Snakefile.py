# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
rule sampling_data_by_tail:
    message:
        '''
        ------------------------------
        sampling raw reads by using tail line of 4 fold 
        ------------------------------
        '''
    input:
        read1 = join(config["samplesDir"], "{sample}"+config["read1Suffix"]),
        read2 = join(config["samplesDir"], "{sample}"+config["read2Suffix"])
    output:
        read1 = join(sampling_data_outdir, "{sample}"+config["read1Suffix"]),
        read2 = join(sampling_data_outdir, "{sample}"+config["read2Suffix"])
    threads:
        5
    benchmark:
        join(sampling_data_outdir, "{sample}", "benchmark.txt")
    log:
        join(sampling_data_outdir, "logfile", "{sample}.log")
    params:
        lines = 40000
    shell:
        '''
        gunzip -c {input.read1} | tail -{params.lines} | gzip > {output.read1} || sleep 0.1 2>{log};
        gunzip -c {input.read2} | tail -{params.lines} | gzip > {output.read2} || sleep 0.1 2>>{log};
        '''