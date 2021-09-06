
rule eggnog:
    input:
        genomeAnno = config['genomeAnno'],
        genomeFasta = onfig['genomeFasta']
    output:
        eggnog_go = ,
        eggnog_kegg = 
    threads: 1
    conda:
        'env/eggnog.yaml'
    log:
       join(),
    shell:
        '''
        emapper.py -i {input.rna} --output {output.dir} -d virNOG --data_dir ../../Database/eggnogdb --translate --cpu {threads} --usemem;
        
        '''