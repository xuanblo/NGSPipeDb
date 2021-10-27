tnt1_outdir = join(config["resultsDir"], "tnt")

rule download_itis_code:
    output:
        code_dir = directory(join(tnt1_outdir, 'software'))
    shell:
        '''
        mkdir -p {output.code_dir};
        cd {output.code_dir};
        git clone git://github.com/Chuan-Jiang/ITIS.git;
        '''

rule tnt1_find_by_itis:
    message:
        '''
        https://github.com/Chuan-Jiang/ITIS

        '''
    input:
        read1 = join(qc_outdir, qc_method, "{sample}", "{sample}.cleanR1.fq.gz"),
        read2 = join(qc_outdir, qc_method, "{sample}", "{sample}.cleanR2.fq.gz"),
        code_dir = join(tnt1_outdir, 'software')
    output:
        filtered_bed = join(tnt1_outdir, 'itis', "{sample}", '{sample}.tnt1.filtered.bed'),
    threads: 30
    params:
        output_dir = directory(join(tnt1_outdir, 'itis', "{sample}"))
    log: join(tnt1_outdir, 'itis', "{sample}", 'run.log.txt')
    shell:
        '''
        #perl {input.code_dir}/ITIS/itis.pl -g {config[genomeFasta]} -t {snake_dir}/metadata/medicago_tnt1.fa -l 300 -N {wildcards.sample} -1 {input.read1} -2 {input.read2} -f {config[genomeAnno]} -c {threads},{threads},{threads} -D {params.output_dir} 1>{log} 2>&1;
        # use -f will may meets $la_no, $ne_no error
        perl {input.code_dir}/ITIS/itis.pl -g {config[genomeFasta]} -t {snake_dir}/metadata/medicago_tnt1.fa -l 300 -N {wildcards.sample} -1 {input.read1} -2 {input.read2} -c {threads},{threads},{threads} -D {params.output_dir} 1>{log} 2>&1;
        '''

rule tnt_merge:
    input:
        output_dir = expand(join(tnt1_outdir, 'itis', '{sample}', '{sample}.tnt1.filtered.bed'), sample=SAMPLES)
    output:
        tnt_merge_ok = touch(join(flag_outdir, 'tnt_merge.ok'))