
rule generate_database:
    output:
        touch(join(config['result_dir'], 'database', 'database.ok'))

rule gtfTosqlite3:
    input:
        genome_gtf = config['genome_gtf']
    output:
        gff_db = join(config['result_dir'], 'database', 'data', 'gtf.sqlite3')
    log: join(config['result_dir'], 'database', 'data', 'gtf.log')
    conda: 'envs/gffutils.yaml'
    shell:
        '''
        #python {source_dirname}/script/gtfTosqlite3.py {input.genome_gtf} {output.gff_db} >{log} 2>&1;
        gffutils-cli create --output {output.gff_db} {input.genome_gtf} >{log} 2>&1
        '''

rule:
    output: touch('okok')
    shell:
        '''
        ncbi-genome-download
        '''

rule anno2sqlite3:
    input:
        annotaion = join(config['result_dir'], 'report', 'anno.tsv')
    output:
        anno_db = join(config['result_dir'], 'database', 'data', 'anno.sqlite3')
    conda: 'envs/simplesqlite.yaml'
    log: join(config['result_dir'], 'database', 'data', 'anno.log')
    shell:
        '''
        pip install SimpleSQLite >{log} 2>&1;
        python {source_dirname}/script/annoTosqlite3.py {input.annotaion} {output.anno_db} >>{log} 2>&1;
        '''


rule blastdb:
    input:
        genome_gtf = config['genome_gtf'],
        genome_fasta = config['genome_fasta'],
        #genome_protein = 'if you have, give a path',
        #genome_cds = 'if you have, give a path',
    output:
        blastdb = directory(join(config['result_dir'], 'database', 'blastdb'))
    conda: 'envs/blast.yaml'
    log: join(config['result_dir'], 'database', 'blastdb', 'makeblasdb.log')
    shell:
        '''
        gffread -g {input.genome_fasta} -G {input.genome_gtf} -x {output}/cds.fa 1>{log} 2>&1;
        gffread -g {input.genome_fasta} -G {input.genome_gtf} -y {output}/protein.fa 1>>{log} 2>&1;
        makeblastdb -in {input.genome_fasta} -dbtype nucl -out {output}/genome 1>>{log} 2>&1;
        makeblastdb -in {output}/protein.fa -dbtype prot -out {output}/protein 1>>{log} 2>&1;
        makeblastdb -in {output}/cds.fa -dbtype nucl -out {output}/cds 1>>{log} 2>&1;
        '''

rule expression_atlas:
    input:
    output:
    shell:
        '''
        '''

rule expression_blast:
    input:
    output:
    shell:
        '''
        '''



