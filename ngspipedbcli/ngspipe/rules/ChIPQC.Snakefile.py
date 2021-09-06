rule macs2_call_peak:    #有control的call peaks
    input:
        treat="bam/10-me_bwa_mm10_sorted_rmdup.bam",    #实验组bam
        control="bam/10-Input_bwa_mm10_sorted_rmdup.bam"    #配对的对照bam
    output:
        "macs2_result/10-me_peaks.broadPeak"
    params:    #params关键词作用，创建局部变量，引用方式 {params.name}
        outdir="macs2_result",
        head_outfile="10-me"
    log:
        "macs2_result/10-me_peaks.log"
    threads:4
    shell:
        "macs2 callpeak -t {input.treat} -c {input.control} -f BAM  -g hs -n {params.head_outfile}  \
         --outdir {params.outdir} \
         --nomodel --extsize 200 \
         -B --broad --broad-cutoff  0.01 "
