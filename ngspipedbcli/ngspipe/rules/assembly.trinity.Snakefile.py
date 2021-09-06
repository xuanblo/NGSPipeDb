# https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html
rule trinity:
    input:
        read1 = join(qc_outdir, "{sample}", "{sample}.cleanR1.fq.gz"),
        read2 = join(qc_outdir, "{sample}", "{sample}.cleanR2.fq.gz"),
    output:
        trinityout = touch(config["resultsDir"] + "/trinity/{sample}/trinity_out.ok")
    params:
        trinityfolder = config["resultsDir"] + "/trinity/{sample}/trinity_out"
    threads: 1
    log:
        trinitysample = config["resultsDir"] + "/trinity/{sample}/run_trinity.log",
    shell:
        '''
        #macos
        #Trinity --seqType fq --max_memory 4G --left {input.read1} --right {input.read2} --output {params.trinityfolder} --CPU {threads} 1>{log} 2>&1;
        ~/Work/Current_work2020-6-21/databasetool/software_mac/trinityrnaseq-v2.12.0/Trinity --seqType fq --max_memory 4G --left {input.read1} --right {input.read2} --output {params.trinityfolder} --CPU {threads} 1>{log} 2>&1;
        '''