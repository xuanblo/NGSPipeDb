rule eggnog:
    input:
        rna = config["resultfolder"] + "/assembly_final/mRNA.fa",
    output:
        dir = directory(config["resultfolder"] + "/eggnog"),
    conda: "env/eggnog.yaml"
    threads: 40
    log: config["resultfolder"] + "/eggnog/run_eggnog.log",
    shell:
        "emapper.py -i {input.rna} --output {output.dir} -d virNOG --data_dir ../../Database/eggnogdb --translate --cpu {threads} --usemem"
