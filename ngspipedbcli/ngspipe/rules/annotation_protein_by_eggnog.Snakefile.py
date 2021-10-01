
rule eggnog_mapper:
    '''
    Document: https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5
    Database:
        download_eggnog_data.py -y --data_dir ./
        tar -zxvf eggnog.db.gz
    '''
    input:
        genomeFasta = config["genomeFasta"],
        genomeAnno = config["genomeAnno"],
    output:
        transcript_fa = join(protein_anno_outdir, 'transcript.fa'),
        eggnog_out = join(protein_anno_outdir, '{}.emapper.annotations'.format('ngspipe_eggnog'))
    params:
        outfile_prefix = 'ngspipe_eggnog',
        tmpdir = join(protein_anno_outdir, '{}.tmp'.format('ngspipe_eggnog'))
    threads: 10
    log: join(protein_anno_outdir, 'run.log.txt')
    shell:
        '''
        if [[ "{config[database_eggnog_dir]}" != "null" && -e "{config[database_eggnog_dir]}" ]]
        then
            gffread -g {input.genomeFasta} {input.genomeAnno} -w {output.transcript_fa} 1>{log} 2>&1;
            mkdir -p {params.tmpdir};
            emapper.py --itype CDS -i {output.transcript_fa} -o {params.outfile_prefix} --output_dir {protein_anno_outdir} -m diamond --cpu {threads} --data_dir {config[database_eggnog_dir]} --override --temp_dir {params.tmpdir} 1>>{log} 2>&1;
        else
            echo "no database"
        fi
        '''

rule eggnog2gokegg:
    '''
    '''
    input:
        eggnog_result = join(protein_anno_outdir, '{}.emapper.annotations'.format('ngspipe_eggnog'))
    output:
        eggnog_GO = join(protein_anno_outdir, 'eggnog.go'),
        eggnog_KO = join(protein_anno_outdir, 'eggnog.ko'),
        eggnog_k = join(protein_anno_outdir, 'eggnog.K'),
    log: join(protein_anno_outdir, 'eggnog2gokegg.log.txt')
    shell:
        '''
        python {snake_dir}/scripts/eggnog2go_kegg.py {config[genomeAnno]} {input.eggnog_result} {output.eggnog_GO} {output.eggnog_KO} {output.eggnog_k} 1>{log} 2>&1;
        '''

rule gokegg_addname:
    '''
    http://yulab-smu.top/clusterProfiler-book/
    '''
    input:
        eggnog_GO = join(protein_anno_outdir, 'eggnog.go'),
        eggnog_KO = join(protein_anno_outdir, 'eggnog.ko'),
        eggnog_k = join(protein_anno_outdir, 'eggnog.K'),
    output:
        eggnog_go_addname = join(protein_anno_outdir, 'eggnog.go.addname.txt'),
        eggnog_kegg_addname = join(protein_anno_outdir, 'eggnog.kegg.addname.txt'),
    log:
        go_addname_log = join(protein_anno_outdir, 'eggnog.go.addname.log.txt'),
        kegg_addname_log = join(protein_anno_outdir, 'eggnog.kegg.addname.log.txt'),
    shell:
        '''
        python {snake_dir}/scripts/add_go_ontology_name.py {config[database_gene_ontology_path]} {input.eggnog_GO} >{output.eggnog_go_addname} 2>{log.go_addname_log};
        Rscript {snake_dir}/scripts/add_kegg_ontology_name.R --k_anno {input.eggnog_k} --output {output.eggnog_kegg_addname} 1>{log.kegg_addname_log} 2>&1;
        '''