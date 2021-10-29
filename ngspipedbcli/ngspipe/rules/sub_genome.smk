sub_genome_dir = join(config["resultsDir"], "sub_genome")

rule sub_genome_by_seqkit:
    message:
        '''
        use part of genome as reference
        '''
    input:
        genomeFasta = config['genomeFasta'],
    output:
        genomeFasta_out = join(sub_genome_dir, 'sub_genome.fa'),
        merge_sub_genome_ok = touch(join(flag_outdir, 'merge_sub_genome.ok'))
    log: join(sub_genome_dir, 'run.log.txt')
    threads: 10
    params:
        cut_method = config['sub_genome']
    shell:
        '''
        python {snake_dir}/scripts/get_sub_genome.py -g {input.genomeFasta} -go {output.genomeFasta_out} -r {params.cut_method} 1>{log} 2>&1;
        '''