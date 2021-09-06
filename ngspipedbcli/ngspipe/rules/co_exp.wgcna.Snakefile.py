rule stat_lncRNA:
    input:
        gtf = config["resultfolder"] + "/cuffcompare/cufcompF.combined.classcode_except_equal.gtf",
    output:
        directory(config["resultfolder"] + "/lncRNA_stat")
    conda: "env/statlncRNA.yaml"
    log: config["resultfolder"] + "/lncRNA_stat/lncRNA_stat.log"
    shell:
        "python script/gff_stat_visualizition_v1.py -i {input.gtf} -o {output} >{log} 2>&1;"
        "python script/filter_lncRNA_stat.py >{output}xxx"