le GRCm38.83.chr19.gtf|awk '$3=="exon"'|awk '{print $1"\t"$4"\t"$5"\texon"NR"\t.\t"$7}' >GRCm38.83.chr19.exon.bed

rule download_sra_data:
    input:
        sralist = config['sralist'],
    output:
        anno = config["resultfolder"] + '/genome_stat/summary.txt',
        fix = config["resultfolder"] + '/genome_stat/stat.txt',
    params:
        outdir = config["resultfolder"] + '/genome_stat'
    log:
        fix = config["resultfolder"] + '/genome_stat/fix.log.txt',
        stat = config["resultfolder"] + '/genome_stat/stat.log.txt',
    conda: "env/s.yaml" # sra-tools
    shell:
        '''
        prefetch --option-file sra.list 
        fastq-dump –X 5 –Z –split-files SRR5907429
        ## 一定要搞清楚你的软件被conda安装在哪
        ## ~/miniconda3/envs/snakemake/etc/asperaweb_id_dsa.openssh
        ## ~/miniconda3/etc/asperaweb_id_dsa.openssh
        ascp -v -k 1 -T -l 200m -i ~/miniconda3/envs/snakemake/etc/asperaweb_id_dsa.openssh dbtest@sra-download.ncbi.nlm.nih.gov:data/sracloud/traces/sra51/SRR/005768/SRR5907429 ./
        for i in `cat sra.list`;do prefetch -a '~/miniconda3/envs/snakemake/bin/ascp|~/miniconda3/envs/snakemake/etc/asperaweb_id_dsa.openssh' $i;done
        prefetch -t fasp -a '/home/zhangxuan/miniconda3/envs/snakemake/bin/ascp|/home/zhangxuan/miniconda3/envs/snakemake/etc/asperaweb_id_dsa.openssh' SRR6656275
        ascp -QT -l 300m -P33001 -i ~/miniconda3/envs/snakemake/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR665/005/SRR6656315/SRR6656315_subreads.fastq.gz .
        wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR665/SRR6656275/SRR6656275_1.fastq.gz .
        le EBI.Study.list|grep PacBio|cut -f10|sed 's/ftp.sra.ebi.ac.uk//'|awk '{print "ascp -QT -P33001 -i ~/miniconda3/envs/snakemake/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:"$1" ."}' >download.sh
        '''


# conda install -y -c hcc aspera-cli
#conda install -y -c bioconda sra-tools
# https://bioinformaticsworkbook.org/dataAcquisition/fileTransfer/sra.html