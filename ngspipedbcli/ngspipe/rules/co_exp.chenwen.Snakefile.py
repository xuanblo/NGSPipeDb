rule coexp:
    input: 
        gene_count = config["resultfolder"] + "/assembly_final/gene.norm.tsv",
    output:
        coexp_network = config["resultfolder"] + "/coexp/co-exp.network",
        stat = config["resultfolder"] + "/coexp/stat.txt",
    log:
        normalized_log = config["resultfolder"] + "/coexp/normalized.log",
        network_construction_log = config["resultfolder"] + "/coexp/network_construction.log",
        network_stat_log = config["resultfolder"] + "/coexp/network_stat.log",
    threads: 40
    shell:
        '''
        perl ../../Software/RNAseq2Coexpnetwork/Final/network_construction.pl --gene_matrix={input.gene_count} --output_network={output.coexp_network} --log2x --cc_method=pcc --cutoff_adjpv=0.05 --cpucore={threads} --outputstat >{log.normalized_log} 2>&1;
        perl ../../Software/RNAseq2Coexpnetwork/Final/network_stat.pl --network={output.coexp_network} --degree_stat={output.stat} >{log.normalized_log} 2>&1;
        '''