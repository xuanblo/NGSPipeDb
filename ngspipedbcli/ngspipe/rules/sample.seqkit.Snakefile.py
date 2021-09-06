rule sample_rawData_by_seqkit:
    input:
        read1 = join(config["sample_dir"], "{sample}_R1.fq.gz"),
        read2 = join(config["sample_dir"], "{sample}_R2.fq.gz")
    output:
        read1 = join(config["result_dir"], "sample.seqkit", "{sample}_R1.fq.gz"),
        read2 = join(config["result_dir"], "sample.seqkit", "{sample}_R2.fq.gz")
    threads: 5
    log: join(config["result_dir"], "sample.seqkit", "{sample}.log"),
    params:
        proportion = 0.1
    conda: "env/pre3.yaml"
    shell:
        '''
        gunzip -c {input.read1} | seqkit sample -p {params.proportion} -j {threads} -o {output.read1} 1>{log} 2>&1;
        gunzip -c {input.read2} | seqkit sample -p {params.proportion} -j {threads} -o {output.read2} 1>>{log} 2>&1;
        '''