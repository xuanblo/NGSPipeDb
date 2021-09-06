# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
rule sampling_data_by_seqkit_proportion:
    message:
        '''
        ------------------------------
        sampling raw reads by using proportion, reads will be selected randomly 
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
    conda:
        "envs/seqkit.yaml"
    params:
        proportion = 0.01
    shell:
        '''
        gunzip -c {input.read1} | seqkit sample -p {params.proportion} -j {threads} -o {output.read1} 1>{log} 2>&1;
        gunzip -c {input.read2} | seqkit sample -p {params.proportion} -j {threads} -o {output.read2} 1>>{log} 2>&1;
        '''

