rule migration_model:
    message:
        '''
        ------------------------------
        8. after model generated, make migration
        ------------------------------
        '''
    input:
        gffdjango_model = join(django_dir, "geneAnno", "models.py"),
        exp_django_model = join(django_dir, "geneExpAtlas", "models.py"),
        copy_ngsdb_file = django_dir,
    output:
        migrationok = touch(join(migration_outdir, "migration.ok")),
    params:
        db_name = "gffDb",
        db_app = 'geneAnno',
    log:
        join(migration_outdir, "change_model_py.log")
    shell:
        '''
        #pip install wooey clustergrammer sklearn pandas==0.25.3 1>{log} 2>&1;
        python {django_dir}/manage.py makemigrations 1>>{log} 2>&1;
        # django makemigrate
        python {django_dir}/manage.py migrate 1>>{log} 2>&1;
        '''
