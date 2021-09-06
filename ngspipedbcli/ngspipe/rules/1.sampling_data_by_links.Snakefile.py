# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
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
        read1 = join(sampling_data_outdir, "{sample}"+config["read1Suffix"]),
        read2 = join(sampling_data_outdir, "{sample}"+config["read2Suffix"])
    threads:
        5
    benchmark:
        join(sampling_data_outdir, "{sample}", "benchmark.txt")
    log:
        join(sampling_data_outdir, "logfile", "{sample}.log")
    run:
        import os
        read1_src = os.path.abspath(input.read1)
        read2_src = os.path.abspath(input.read2)
        os.symlink(read1_src, output.read1)
        os.symlink(read2_src, output.read2)