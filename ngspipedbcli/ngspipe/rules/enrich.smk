
gokegg_enrich_outdir = join(config["resultsDir"], "gokegg_enrich")
enrich_method = config['enrich_method']

rule diff_gene_enrich_by_clusterprofiler:
    input:
        diff_outputok = join(flag_outdir, 'differential_expression.ok'),
        eggnog_go_addname = join(protein_anno_outdir, anno_method, 'eggnog.go.addname.txt'),
        eggnog_kegg_addname = join(protein_anno_outdir, anno_method, 'eggnog.kegg.addname.txt'),
    output:
         diff_gene_enrich_ok = touch(join(gokegg_enrich_outdir, 'clusterprofiler', 'diff_gene_enrich_ok')),
    log: join(gokegg_enrich_outdir, 'clusterprofiler', 'run.log.txt')
    params:
        diff_gene_enrich_outputdir = join(gokegg_enrich_outdir, 'clusterprofiler')
    shell:
        '''
        mkdir -p {params.diff_gene_enrich_outputdir};
        for i in `ls {diff_outdir}/{diff_method}/*.padj.csv|sed 's#.*/##'|sed 's#.csv##'`;do Rscript {snake_dir}/scripts/run_diff_gene_go_enrich.R --difflist {diff_outdir}/{diff_method}/$i.csv --term_anno {input.eggnog_go_addname} --output {params.diff_gene_enrich_outputdir}/$i.go_enrich.csv;done 1>{log} 2>&1;
        for i in `ls {diff_outdir}/{diff_method}/*.padj.csv|sed 's#.*/##'|sed 's#.csv##'`;do Rscript {snake_dir}/scripts/run_diff_gene_kegg_enrich.R --difflist {diff_outdir}/{diff_method}/$i.csv --term_anno {input.eggnog_kegg_addname} --output {params.diff_gene_enrich_outputdir}/$i.kegg_enrich.csv;done 1>>{log} 2>&1;
        '''

rule diff_gene_enrich_by_touch_empty:
    input:
        diff_outputok = join(flag_outdir, 'differential_expression.ok'),
        eggnog_go_addname = join(protein_anno_outdir, anno_method, 'eggnog.go.addname.txt'),
        eggnog_kegg_addname = join(protein_anno_outdir, anno_method, 'eggnog.kegg.addname.txt'),
    output:
        diff_gene_enrich_ok = touch(join(gokegg_enrich_outdir, 'touch_empty', 'diff_gene_enrich_ok')),

rule enrich:
    message:
        '''
        differential_expression ok
        '''
    input:
        diff_gene_enrich_ok = join(gokegg_enrich_outdir, enrich_method, 'diff_gene_enrich_ok'),
        protein_annotation_ok = join(flag_outdir, 'protein_annotation.ok')
    output:
        enrich_ok = touch(join(flag_outdir, 'enrich.ok'))
