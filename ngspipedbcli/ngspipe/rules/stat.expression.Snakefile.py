le GRCm38.83.chr19.gtf|awk '$3=="exon"'|awk '{print $1"\t"$4"\t"$5"\texon"NR"\t.\t"$7}' >GRCm38.83.chr19.exon.bed

rule exp_stat:
    input:
        gene_exp_matrix = config["resultfolder"] + "/assembly_final/gene_fpkm_all_samples.tsv",
        lncRNA_id_file = config["resultfolder"] + "/assembly_final/lncRNA.fa",
        mRNA_id_file = config["resultfolder"] + "/assembly_final/mRNA.fa",
    output:
        fpkm_distribution = config["resultfolder"] + "/exp_stat/fpkm_distribution.pdf"
    shell:
        '''
        python script/exp_stat.py {input.gene_exp_matrix} {input.lncRNA_id_file} {input.mRNA_id_file} {output.fpkm_distribution}
        '''