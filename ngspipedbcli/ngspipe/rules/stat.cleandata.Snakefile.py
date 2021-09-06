rule clean_data_stat:
    input:
        fwd = config["datafolder"] + "/{sample}_1.clean.fq.gz",
        rev = config["datafolder"] + "/{sample}_2.clean.fq.gz"
    output:
        fwd = config["resultfolder"] + "/clean_data_stat/{sample}/{sample}_1.seqkit.summary",
        rev = config["resultfolder"] + "/clean_data_stat/{sample}/{sample}_2.seqkit.summary",
        fwd_gc_qual = config["resultfolder"] + "/clean_data_stat/{sample}/{sample}_1.seqkit.gc_qual",
        rev_gc_qual = config["resultfolder"] + "/clean_data_stat/{sample}/{sample}_2.seqkit.gc_qual",
    threads: 1
    shell:
        '''
        seqkit stat -j {threads} -a {input.fwd} >{output.fwd};
        seqkit stat -j {threads} -a {input.rev} >{output.rev};
        seqkit fx2tab -n -g -q -j {threads} {input.fwd}|cut -f4-|awk '{{gc=gc+$1;qual=qual+$2;}}END{{print gc/NR"\\t"qual/NR}}' >{output.fwd_gc_qual};
        seqkit fx2tab -n -g -q -j {threads} {input.rev}|cut -f4-|awk '{{gc=gc+$1;qual=qual+$2;}}END{{print gc/NR"\\t"qual/NR}}' >{output.rev_gc_qual};
        '''



# ----------------------------------------------------------------------------
# refer_library_by_infer_experiment_RseQC
multiQC_ourdir = join(config["result_dir"], "qc.multiQC")
# ----------------------------------------------------------------------------
rule QC_by_multiQC:
    input:
        expand(join(trim_galore_ourdir, "{sample}", "{sample}_val_1_fastqc.zip"), sample=SAMPLES),
        expand(join(trim_galore_ourdir, "{sample}", "{sample}_val_2_fastqc.zip"), sample=SAMPLES)
    output: touch(join(multiQC_ourdir, "multiqc.ok"))
    conda: "env/pre3.yaml"
    log: join(multiQC_ourdir, "multiqc.log")
    shell:
        '''
        multiqc -f -o {multiQC_ourdir} -n multiqc {trim_galore_ourdir} 1>{log} 2>&1;
        '''