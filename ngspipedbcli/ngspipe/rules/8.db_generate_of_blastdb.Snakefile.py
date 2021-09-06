genomedb_outdir = join(blastdb_outdir, 'nucl_genomedb')
transcriptdb_outdir = join(blastdb_outdir, 'nucl_transcriptdb')
proteindb_outdir = join(blastdb_outdir, 'prot_proteindb')

rule merge_blastdb:
    message:
        '''
        ------------------------------
        8. merge all blast database
        ------------------------------
        '''
    input:
        genomeFastaDb = join(genomedb_outdir, "genomefasta_blastdb.ok")
    output:
        merged_db = touch(join(blastdb_outdir, "merged_blastdb.ok"))

rule makeblastdb_genomefasta_by_blast:
    message:
        '''
        ------------------------------
        8. generate genome blast database
        ------------------------------
        '''
    input:
        genomeFasta = config['genomeFasta']
    output:
        genomeFastaDb = touch(join(genomedb_outdir, "genomefasta_blastdb.ok"))
    log:
        join(genomedb_outdir, "run.log")
    benchmark:
        join(genomedb_outdir, "benchmark.txt")
    conda:
        '../envs/blast.yaml'
    params:
        genomedb_prefix = 'genome'
    shell:
        '''
        makeblastdb -in {input.genomeFasta} -dbtype nucl -out {genomedb_outdir}/{params.genomedb_prefix} 1>{log} 2>&1;
        '''


rule makeblastdb_transcriptfasta_by_blast:
    message:
        '''
        ------------------------------
        8. generate transcript blast database
        ------------------------------
        '''
    input:
        genomeFasta = config['genomeFasta'],
        genomeAnno = config['genomeAnno'],
    output:
        transcriptFastaDb = touch(join(transcriptdb_outdir, "nucl", "transcriptfasta_blastdb.ok"))
    log:
        join(transcriptdb_outdir, "run.log")
    benchmark:
        join(transcriptdb_outdir, "benchmark.txt")
    conda:
        '../envs/blast.yaml'
    params:
        genomedb_prefix = 'genome'
    shell:
        '''
        gffread {input.genomeAnno} -g {input.genomeFasta} -w {output.transcriptdb_outdir}/nucl/transcript.fa
        makeblastdb -in {output.transcriptdb_outdir}/nucl/transcript.fa -dbtype nucl -o {blastdb_outdir}/{params.genomedb_prefix} 1>{log} 2>&1;
        '''
