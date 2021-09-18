rule interproscan:
    '''
    Document: https://interproscan-docs.readthedocs.io/en/latest/
    couldn't install by conda
    '''
    input:
        genomeFasta = config["genomeFasta"],
        genomeAnno = config["genomeAnno"],
    output:
        transcript_fa = join(protein_anno_outdir, 'transcript.fa'),
        interproscan_outdir = directory(join(protein_anno_outdir, 'interproscan'))
    log: join(protein_anno_outdir, 'run.log.txt')
    threads: 10
    shell:
        '''
        wget http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.52-86.0/interproscan-5.52-86.0-64-bit.tar.gz && tar -zxvf interproscan-5.52-86.0-64-bit.tar.gz;
        interproscan.sh -cpu {threads} --formats TSV,html --goterms --highmem --pathways --tempdir {output.outdir}/tmpfile -iprlookup -i {input.genom
e_pep} -d {output.outdir} 1>{log} 2>&1;
        '''