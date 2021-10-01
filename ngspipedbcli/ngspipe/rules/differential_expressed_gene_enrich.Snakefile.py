rule diff_gene_enrich:
    input:
        diff_outputok = join(diff_outdir, "diff.ok"),
        eggnog_go_addname = join(protein_anno_outdir, 'eggnog.go.addname.txt'),
        eggnog_kegg_addname = join(protein_anno_outdir, 'eggnog.kegg.addname.txt'),
    output:
         diff_gene_enrich_outputdir = directory(join(gokegg_enrich_outdir)),
    log: join(gokegg_enrich_outdir, 'run.log.txt')
    shell:
        '''
        mkdir -p {gokegg_enrich_outdir}
        for i in `ls {diff_outdir}/*.padj.csv|sed 's#.*/##'|sed 's#.csv##'`;do Rscript {snake_dir}/scripts/run_diff_gene_go_enrich.R --difflist {diff_outdir}/$i.csv --term_anno {input.eggnog_go_addname} --output {output.diff_gene_enrich_outputdir}/$i.go_enrich.csv;done 1>{log} 2>&1;
        for i in `ls {diff_outdir}/*.padj.csv|sed 's#.*/##'|sed 's#.csv##'`;do Rscript {snake_dir}/scripts/run_diff_gene_kegg_enrich.R --difflist {diff_outdir}/$i.csv --term_anno {input.eggnog_kegg_addname} --output {output.diff_gene_enrich_outputdir}/$i.kegg_enrich.csv;done 1>>{log} 2>&1;
        '''
