rule cutadapt:    #cutadapt关键词，定义cutdapt的流程，主要包括input，output，shell三个部分
    input:
        raw_R1 = join(config['samplesDir'], "{sample}"+config["read1Suffix"]),    #rawdata的R1.fastq
        raw_R2 = join(config['samplesDir'], "{sample}"+config["read2Suffix"])    #rawdata的R2.fastq
    output:
        clean_R1 = join(config['resultsDir'], "chipseq", "clean_fastq", "{sample}"+config["read1Suffix"]),    #cutadapt输出的R1.fastq
        clean_R2 = join(config['resultsDir'], "chipseq", "clean_fastq", "{sample}"+config["read2Suffix"]),    #cutadapt输出的R2.fastq
    log:
        "clean_fastq/{sample}_cutadapt.log"
    threads: 4    #所需的核心数，与-j保持一致
    shell:    #cutadapt的代码
        "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -u 5 -u -0 -U 5 -U -0 -m 30 -j 4 -o {output.clean_R1} -p {output.clean_R2} {input.raw_R1} {input.raw_R2} 1>{log} 2>&1"
