rule download_sra_by_grabseqs:
    input:
    output:
    conda: 'sra-tools.yaml'
    shell:
        '''
        pip install grabseqs;
        grabseqs sra -o 黑色素瘤 SRS8434327
        '''

rule download_sra_by_sra-tools:
    input:
    output:
    shell:
        '''
        prefetch --option-file SraAccList.txt
        prefetch -O 黑色素瘤/ SRS8434327 -p
        prefetch -O 黑色素瘤/ SRS8434327 -p --ascp-path "~/miniconda3/envs/sradownload/bin/ascp|~/miniconda3/envs/sradownload/etc/asperaweb_id_dsa.openssh"
        '''
