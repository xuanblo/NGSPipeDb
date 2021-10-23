network_method = config['network_method']
rule network_analysis_by_gcen:
    input:
        normalized_counts = join(diff_outdir, diff_method, 'normalized.counts.csv'),
    output:
        network = join(coexp_outdir, 'gcen', 'co_exp.network.tsv')
    threads: 20
    log: join(coexp_outdir, 'gcen', 'run.log.txt')
    shell:
        '''
        less {input.normalized_counts} |tr ',' '\\t' >{coexp_outdir}/exp.tsv
        network_build -i {coexp_outdir}/exp.tsv -o {output.network} -t {threads} -p 1 -c 0 1>{log} 2>&1;
        '''

rule network:
    input:
        network = join(coexp_outdir, network_method, 'co_exp.network.tsv'),
        enrich_ok = join(flag_outdir, 'enrich.ok')
    output:
        network_ok = touch(join(flag_outdir, 'network.ok'))