rule timeSerise:
    input:
        gene_exp_matrix = config["resultfolder"] + "/assembly_final/gene_fpkm_all_samples.tsv",
    output:
        gene_clusters = config["resultfolder"] + "/timeSerise/result.txt",
    params:
        ourdir = config["resultfolder"] + "/timeSerise"
    conda: "env/clust.yaml"
    shell:
        '''
        # IB1->IB2->IB3->FFB->FF
        python script/sub_exp_matrix.py IB1:IB2:IB3:FFB:FF {input.gene_exp_matrix} {params.ourdir}/IB1_IB2_IB3_FFB_FF.matrix;
        python script/exp_matrix2replicate.py {params.ourdir}/IB1_IB2_IB3_FFB_FF.matrix >{params.ourdir}/IB1_IB2_IB3_FFB_FF.txt;
        # IB1->IB2->IB3->MFB->MF
        python script/sub_exp_matrix.py IB1:IB2:IB3:MFB:MF {input.gene_exp_matrix} {params.ourdir}/IB1_IB2_IB3_MFB_MF.matrix;
        python script/exp_matrix2replicate.py {params.ourdir}/IB1_IB2_IB3_MFB_MF.matrix >{params.ourdir}/IB1_IB2_IB3_MFB_MF.txt;
        # cluster
        clust {params.ourdir}/IB1_IB2_IB3_FFB_FF.matrix -r {params.ourdir}/IB1_IB2_IB3_FFB_FF.txt -o {params.ourdir}/clust_result_IB1_IB2_IB3_FB_FF;
        clust {params.ourdir}/IB1_IB2_IB3_MFB_MF.matrix -r {params.ourdir}/IB1_IB2_IB3_MFB_MF.txt -o {params.ourdir}/clust_result_IB1_IB2_IB3_MFB_MF;
        '''

rule tissueEnrich:
    input:
        gene_exp_matrix = config["resultfolder"] + "/assembly_final/gene_fpkm_all_samples.tsv",
    output: 
        tissueSpecified_gene = config["resultfolder"] + "/tissueEnrich/ts.tsv",
    conda: "env/tissueEnrich.yaml"
    shell:
        '''
        Rscript script/TissueEnrich.R {input.gene_exp_matrix} {output.tissueSpecified_gene};
        '''
