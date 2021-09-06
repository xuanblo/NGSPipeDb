le GRCm38.83.chr19.gtf|awk '$3=="exon"'|awk '{print $1"\t"$4"\t"$5"\texon"NR"\t.\t"$7}' >GRCm38.83.chr19.exon.bed

rule genome_stat:
    input:
        genome = config['genome'],
        gff = config['gff']
    output:
        anno = config["resultfolder"] + '/genome_stat/summary.txt',
        fix = config["resultfolder"] + '/genome_stat/stat.txt',
    params:
        outdir = config["resultfolder"] + '/genome_stat'
    log:
        fix = config["resultfolder"] + '/genome_stat/fix.log.txt',
        stat = config["resultfolder"] + '/genome_stat/stat.log.txt',
        statistics = config["resultfolder"] + '/genome_stat/statistics.log.txt',
        dist = config["resultfolder"] + '/genome_stat/dist.log.txt',
    conda: "env/gaas.yaml"
    shell:
        '''
            # fix or repair gff file
    agat_sp_gxf_to_gff3.pl -g xingyouteng.gene.final.changeid.gff -o Pvo.gff
    agat_sp_add_start_and_stop.pl --gff Pvo.gff --fasta xingyouteng.FINAL.fasta --out Pvo1.gff
    agat_sp_add_introns.pl --gff Pvo1.gff --out Pvo2.gff
    perl gff3sort/gff3sort.pl --chr_order natural Pvo2.gff >Pvo3.gff
            perl {params.gaasBin}/gxf_to_gff3.pl -g {input.gff} -o {params.outdir}/fixed.gff 1>{log.fix} 2>&1;
            perl {params.gaasBin}/gff3_sp_add_introns.pl --gff={params.outdir}/fixed.gff --out={params.outdir}/fixed.intron 1>{log.stat} 2>&1;
            # statistic
            perl {params.gaasBin}/gff3_sp_statistics.pl --gff={params.outdir}/fixed.gff -o {params.outdir}/fixed.statistics.txt 1>{log.statistics} 2>&1;
            python script/genome_stat2csv.py {params.outdir}/fixed.statistics.txt {params.outdir}/fixed.statistics.csv 
            # distribution
            python script/gff_stat_visualizition_v2.py -i {params.outdir}/fixed.intron.gff -g {input.genome} -o {params.outdir} 1>{log.dist} 2>&1;
            # le ../../Database/JcGenome/GCF_000696525.1_JatCur_1.0/GCF_000696525.1_JatCur_1.0_genomic.gff|grep -E 'mRNA|protein_coding' >../../Database/JcGenome/GCF_000696525.1_JatCur_1.0/GCF_000696525.1_JatCur_1.0_mRNA.gff
            # gffread ../../Database/JcGenome/GCF_000696525.1_JatCur_1.0/GCF_000696525.1_JatCur_1.0_mRNA.gff -T -o ../../Database/JcGenome/GCF_000696525.1_JatCur_1.0/GCF_000696525.1_JatCur_1.0_mRNA.gtf
        '''
