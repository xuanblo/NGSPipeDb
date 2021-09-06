rule sqlite3_auto_by_expression_data:
    message:
        '''
        ------------------------------
        8. generate expression data sqlite3
        ------------------------------
        '''
    input:
        gene_fpkm = join(quantify_outdir, "gene_fpkm_all_samples.tsv")
    output:
        exp_sqlite3 = join(exp_db_outdir, "exp.sqlite3")
    log:
        join(exp_db_outdir, "run.log")
    benchmark:
        join(exp_db_outdir, "benchmark.txt")
    shell:
        '''
        #pip install simplesqlite 1>{log} 2>&1;
        python {snake_dir}/scripts/exp2sqlite3.py {input.gene_fpkm} {output.exp_sqlite3} 1>>{log} 2>&1;
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
        #pip install SimpleSQLite >{log} 2>&1;
        python {snake_dir}/script/annoTosqlite3.py {input.annotaion} {output.anno_db} >>{log} 2>&1;
        '''