rule network_analysis_by_gcen:
    input:
        normalized_counts = join(diff_outdir, 'normalized.counts.csv'),
    output:
        network = join(coexp_outdir, 'co_exp.network.tsv')
    threads: 20
    log: join(coexp_outdir, 'run.log.txt')
    shell:
        '''
        less {input.normalized_counts} |tr ',' '\\t' >{coexp_outdir}/exp.tsv
        network_build -i {coexp_outdir}/exp.tsv -o {output.network} -t {threads} -p 1 -c 0 1>{log} 2>&1;
        '''