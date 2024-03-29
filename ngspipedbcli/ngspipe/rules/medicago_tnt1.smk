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
        code_dir = join(tnt1_outdir, 'software'),
        genomeFasta_out = join(sub_genome_dir, 'sub_genome.fa'),
        exogenous_seq = config['exogenous_seq_path']
    output:
        #filtered_bed = join(tnt1_outdir, 'itis', "{sample}", '{sample}.tnt1.filtered.bed'),
        #empty = join(tnt1_outdir, 'itis', "{sample}", '{sample}.tnt1.empty'),
        runned = touch(join(tnt1_outdir, 'itis', "{sample}", '{sample}.runned')),
    threads: 30
    params:
        output_dir = directory(join(tnt1_outdir, 'itis', "{sample}"))
    log: join(tnt1_outdir, 'itis', "{sample}", 'run.log.txt')
    shell:
        '''
        # use -f will may meets $la_no, $ne_no error
        perl {input.code_dir}/ITIS/itis.pl -g {input.genomeFasta_out} -t {input.exogenous_seq} -l 300 -N {wildcards.sample} -1 {input.read1} -2 {input.read2} -f {config[genomeAnno]} -c {threads},{threads},{threads} -D {params.output_dir} 1>{log} 2>&1;
        '''

rule tnt_merge:
    input:
        runned = expand(join(tnt1_outdir, 'itis', '{sample}', '{sample}.runned'), sample=SAMPLES)
    output:
        tnt_merge_ok = touch(join(flag_outdir, 'tnt_merge.ok'))