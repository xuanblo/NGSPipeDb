# ----------------------------------------------------------------------------
# variation annotation

# ----------------------------------------------------------------------------
rule genome_build_by_snpEff:
    message:
        '''
        ------------------------------
        snpEff
        ------------------------------
        '''
    input:
        genomeFasta = config["genomeFasta"],
        genomeAnno = config['genomeAnno'],
    output:
        configfile = join(annotation_outdir, "snpeff_genomebuild", "snpEff.config"),
        genomeok = touch(join(annotation_outdir, "snpeff_genomebuild", "build.ok")),
    params:
        dbdir = join(annotation_outdir, "snpeff_genomebuild", "data"),
        dbdir_relative2configfile = 'data',
        dbname = dbname,
        annoFormat = annoFormat # gff3, gff2, gtf22
    benchmark:
        join(annotation_outdir, "genomebuild", "benchmark.txt")
    log:
        join(annotation_outdir, "genomebuild", "build.log")
    threads:
        5
    shell:
        '''
        # 1. Edit the config file to create the new genome:
        echo "# Third party databases" >{output.configfile} 2>{log};
        echo "{params.dbname}.genome : resequencing local genome" >>{output.configfile} 2>>{log};
        echo "data.dir = {params.dbdir_relative2configfile}" >>{output.configfile} 2>>{log};
        # 2. cp genome files
        mkdir -p {params.dbdir}/genomes 1>>{log} 2>&1;
        mkdir -p {params.dbdir}/{params.dbname} 1>>{log} 2>&1;
        cp {input.genomeFasta} {params.dbdir}/genomes/{params.dbname}.fa 1>>{log} 2>&1;
        cp {input.genomeAnno} {params.dbdir}/{params.dbname}/genes.gtf 1>>{log} 2>&1;
        # 3. build
        snpEff build -{params.annoFormat} -c {output.configfile} {params.dbname} 1>>{log} 2>&1;
        '''

rule variation_annotation_by_snpEff:
    message:
        '''
        ------------------------------
        snpEff
        ------------------------------
        '''
    input:
        genomeFasta = config["genomeFasta"],
        vcf = join(variation_outdir, "{sample}", "{sample}.vcf"),
        genomeok = join(annotation_outdir, "snpeff_genomebuild", "build.ok"),
        configfile = join(annotation_outdir, "snpeff_genomebuild", "snpEff.config"),
    output:
        eff = join(annotation_outdir, "{sample}", "{sample}.snpeff.vcf"),
    params:
        dbname = dbname,
        summary = join(annotation_outdir, "{sample}", "{sample}.snpeff.summary.html"),
    benchmark:
        join(annotation_outdir, "{sample}", "benchmark.txt")
    log:
        join(annotation_outdir, "{sample}", "effect.log")
    threads:
        5
    shell:
        '''
        snpEff eff -c {input.configfile} -s {params.summary} -o vcf {params.dbname} {input.vcf} >{output.eff} 2>{log};
        '''


