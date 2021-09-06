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
        dna2rna_path = join(django_dir, 'wooey', 'wooey_scripts', 'dna2rna.py'),
        dna2rna_name = 'dna2rna_convert',
    shell:
        '''
        # remove scripts in wooey database
        python {snake_dir}/scripts/remove_script_in_wooey_db.py {django_dir} 1>{log} 2>&1
        # add scripts in current machine
        python {django_dir}/manage.py addscript {params.dna2rna_path} --name {params.dna2rna_name} 1>>{log} 2>&1;
        '''