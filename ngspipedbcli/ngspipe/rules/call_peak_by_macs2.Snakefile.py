rule macs2_call_peak:
    input:
        rmdup_bam = expand(join(remove_duplication_outdir, "{sample}", "{sample}.sorted_rmdup.bam"), sample=SAMPLES),
    output:
        callok = touch(join(call_peak_outdir, "callpeak.ok"))
    params:
        outdir = join(call_peak_outdir, "peaks.broadPeak"),
        out_prefix = "10-me"
    log:
        join(call_peak_outdir, "run.log")
    threads:4
    shell:
        '''
        macs2 callpeak -t {input.rmdup_bam[1]} -c {input.rmdup_bam[0]} -f BAMPE -g hs -n {params.out_prefix} --outdir {params.outdir} --keep-dup auto --bdg 1>{log} 2>&1;
        '''


