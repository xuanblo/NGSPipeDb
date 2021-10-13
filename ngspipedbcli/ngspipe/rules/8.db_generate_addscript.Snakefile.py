rule add_script2wooey:
    message:
        '''
        ------------------------------
        8. add_script2wooey
        ------------------------------
        '''
    input:
        makemigration                   = join(migration_outdir, "migration.ok"),
    output:
        addscript                       = touch(join(addscript_outdir, "addscript.ok")),
    log:
        join(working_dir, addscript_outdir, "run.log")
    benchmark:
        join(working_dir, addscript_outdir, "benchmark.txt")
    params:
        # script 1
        dna2rna_path = join(django_dir, 'wooey', 'wooey_scripts', 'dna2rna.py'),
        dna2rna_name = 'dna2rna_convert',
        # script 2
        goenrich_path = join(django_dir, 'wooey', 'wooey_scripts', 'goenrich.py'),
        goenrich_name = 'goenrich',
        # script 3
        keggenrich_path = join(django_dir, 'wooey', 'wooey_scripts', 'keggenrich.py'),
        keggenrich_name = 'keggenrich',
        # script 4
        getsequence_path = join(django_dir, 'wooey', 'wooey_scripts', 'getsequence.py'),
        getsequence_name = 'getsequence',
    shell:
        '''
        # remove scripts in wooey database
        python {snake_dir}/scripts/remove_script_in_wooey_db.py {django_dir} 1>{log} 2>&1
        # add scripts in current machine
        python {django_dir}/manage.py addscript {params.goenrich_path} --name {params.goenrich_name} 1>>{log} 2>&1;
        python {django_dir}/manage.py addscript {params.dna2rna_path} --name {params.dna2rna_name} 1>>{log} 2>&1;
        python {django_dir}/manage.py addscript {params.keggenrich_path} --name {params.keggenrich_name} 1>>{log} 2>&1;
        python {django_dir}/manage.py addscript {params.getsequence_path} --name {params.getsequence_name} 1>>{log} 2>&1;
        '''