# ----------------------------------------------------------------------------
#
stringtieMerge_outdir = join(config["result_dir"], "align.hisat2_stringtie", "step4_merge_stringtieResult_by_stringtieMerge")
feelnc_outdir = join(config["result_dir"], "coding_potential.all", "feelnc")
stringtieMerge_outdir_cuffPrefix = "cufcompF"
# ----------------------------------------------------------------------------
rule feelnc:
    input: 
        refgtf = config["genome_gff"],
        candidaregtf = join(stringtieMerge_outdir, "cufcompF.combined.gtf"),
        genome = config["genome_fasta"]
    output:
        classify = join(feelnc_outdir, "lncRNA_classes.txt"),
    threads: 20
    log:
        filter_log = config["resultfolder"] + "/feelnc/feelnc.log"
    conda: "env/feelnc.yaml"
    shell:
        '''
        python script/repair_gff.py -i {input.refgtf} -o {params.outdir}/repaired.gff;
        awk '$3==\"exon\"||$3==\"mRNA\"' {params.outdir}/repaired.gff|grep -P 'gbkey=mRNA' >{params.outdir}/repaired2.gff;
        gffread {params.outdir}/repaired2.gff -F -T -o {params.outdir}/repaired3.gtf;
        less {params.outdir}/repaired3.gtf|sed 's/transcript_id/;transcript_id/'|awk -F';' '{{print $1$3\";\"$2}}'|sed 's/ gene_id/gene_id/'|sed 's/;/; /' >{params.outdir}/mRNA.gtf;
        # filter
        FEELnc_filter.pl -i {input.candidaregtf} --mRNAfile {params.outdir}/mRNA.gtf --monoex=-1 --size=200 -p {threads} -o {log.filter_log} >{params.outdir}/candidate_lncRNA.gtf 2>{log.filter_log}.filter;
        # coding potential
        FEELnc_codpot.pl -i {params.outdir}/candidate_lncRNA.gtf -a {params.outdir}/mRNA.gtf -g {input.genome} --mode=shuffle --outdir {params.outdir}/codpot 2>{log.filter_log}.codpot;
        # class
        FEELnc_classifier.pl -i {params.outdir}/codpot/candidate_lncRNA.gtf.lncRNA.gtf -a {params.outdir}/mRNA.gtf -l {log.filter_log}.classifier > {params.outdir}/lncRNA_classes.txt ;
        gffread {params.outdir}/codpot/candidate_lncRNA.gtf.lncRNA.gtf -g {input.genome} -w {params.outdir}/candidate_lncRNA.fa;
        FEELnc_pipeline.sh --candidate=result/cuffcompare/cufcompF.combined.gtf --reference=result/feelnc/repaired.gtf --genome=../../Database/JcGenome/GCF_000696525.1_JatCur_1.0/GCF_000696525.1_JatCur_1.0_genomic.fna --outname=FEELnc_pipeline --outdir=result/feelnc >{log} 2>&1;
        '''

rule filtermRNA:
    input:
        gtf = config["resultfolder"] + "/feelnc/candidate_lncRNA.fa",
        genome = config["genome"]
    output:
        gtf = config["resultfolder"] + "/filtermRNA/cufcompF.combined.uix.gtf",
        fasta = config["resultfolder"] + "/filtermRNA/cufcompF.combined.uix.fa"
    conda: "env/gffread.yaml"
    shell:
        "script/filter_mRNA_step0.sh {input.gtf} >{output.gtf};"
        "gffread {output.gtf} -g {input.genome} -w {output.fasta};"
        
# ----------------------------------------------------------------------------
#
cpc2_outdir = join(config["result_dir"], "coding_potential.all", "cpc2")
# ----------------------------------------------------------------------------
rule cpc2:
    input:
        join(stringtieMerge_outdir, "cufcompF.combined.classcode_except_equal.fa")
    output:
        join(cpc2_outdir, "cpc2.out.txt")
    conda:
        "env/pre2.yaml"
    log:
        join(cpc2_outdir, "cpc2.log")
    shell:
        '''
        python ../../Software/CPC2-beta/bin/CPC2.py -i {input} -o {output} >{log} 2>&1;
        '''

# ----------------------------------------------------------------------------
#
cnci_outdir = join(config["result_dir"], "coding_potential.all", "cnci")
# ----------------------------------------------------------------------------
rule cnci:
    input: config["resultfolder"] + "/cuffcompare/cufcompF.combined.classcode_except_equal.fa"
    output: config["resultfolder"] + "/codingPotential/CNCI.index"
    conda: "env/pre2.yaml"
    params:
        outdir = config["resultfolder"] + "/codingPotential"
    threads: 40
    log: config["resultfolder"] + "/codingPotential/cnci.log"
    shell:
        '''
        cd ../../Software;
        python CNCI_package/CNCI.py -f ../Analysis/RunAll/{input} -o ../Analysis/RunAll/{params.outdir} -m pl -p {threads} >../Analysis/RunAll/{log} 2>&1;
        '''

# plek change sys.exit(1) to sys.exit(0)
# ----------------------------------------------------------------------------
#
PLEK_outdir = join(config["result_dir"], "coding_potential.all", "PLEK")
# ----------------------------------------------------------------------------
rule plek:
    input: config["resultfolder"] + "/feelnc/candidate_lncRNA.fa"
    output: config["resultfolder"] + "/codingPotential/PLEK_out"
    threads: 40
    log: config["resultfolder"] + "/codingPotential/PLEK.log"
    conda: "env/plek.yaml"
    shell:
        '''
        python ../../Software/PLEK.1.2/PLEK.py -fasta {input} -out {output} -thread {threads} >{log} 2>&1; # line 407 sys.exit(0)
        '''

# https://pypi.org/project/gff3tool/
# ----------------------------------------------------------------------------
#
pfam_outdir = join(config["result_dir"], "coding_potential.all", "pfam")
# ----------------------------------------------------------------------------
rule pfam:
    input: config["resultfolder"] + "/feelnc/candidate_lncRNA.fa"
    output: config["resultfolder"] + "/pfam/pfamscan.txt"
    params:
        pfamDB = "../../Database/pfam"
    threads: 40
    conda: "env/pfam.yaml"
    log: config["resultfolder"] + "/pfam/run_pfam.log"
    shell:
        "pfam_scan.pl -fasta {input} -dir {params.pfamDB} -out {output} -cpu {threads} >{log} 2>&1"

# ----------------------------------------------------------------------------
#
mergeCodingPotential_outdir = join(config["result_dir"], "coding_potential.all")
# ----------------------------------------------------------------------------
rule mergeCodingPotential:
    input:
        cpc2 = config["resultfolder"] + "/codingPotential/cpc2.out.txt",
        cnci = config["resultfolder"] + "/codingPotential/CNCI.index",
        plek = config["resultfolder"] + "/codingPotential/PLEK_out",
        pfam = config["resultfolder"] + "/pfam/pfamscan.txt",
        mRNA_gtf = config["resultfolder"] + "/feelnc/mRNA.gtf",
    output:
        tlist = config["resultfolder"] + "/codingPotential/noncoding_union.tlist.txt",
    conda: "env/matplotlib.yaml"
    log: config["resultfolder"] + "/codingPotential/mergeCodingPotential.log"
    params:
        outdir = config["resultfolder"] + "/codingPotential",
        feelncdir = config["resultfolder"] + "/feelnc",
    shell:
        '''
        python script/mergeCodingPotential.py -c {input.cpc2} -n {input.cnci} -p {input.plek} --pfam {input.pfam} -o {params.outdir} >{log} 2>&1;
        for i in `cat {params.outdir}/noncoding_union.tlist.txt <(cut -f3 {params.feelncdir}/lncRNA_classes.txt|grep -v 'lncRNA_transcript'|sort -u)|sort|uniq -d`;do grep $i {params.feelncdir}/candidate_lncRNA.gtf;done >{params.outdir}/noncoding_union.gtf;
        cat {params.outdir}/noncoding_union.gtf {input.mRNA_gtf} >{params.outdir}/transcript_candidate.gtf;
        '''

rule mergeCodingPotential:
    input:
        cpc2 = config["resultfolder"] + "/codingPotential/cpc2.out.txt",
        cnci = config["resultfolder"] + "/codingPotential/CNCI.index",
        plek = config["resultfolder"] + "/codingPotential/PLEK_out",
        pfam = config["resultfolder"] + "/pfam/pfamscan.txt",
        candidate_gtf = config["resultfolder"] + "/cuffcompare/cufcompF.combined.classcode_except_equal.gtf",
        equal_gtf = config["resultfolder"] + "/cuffcompare/cufcompF.combined.classcode_equal.gtf",
    output:
        tlist = config["resultfolder"] + "/codingPotential/noncoding_union.tlist.txt",
    conda: "env/matplotlib.yaml"
    log: config["resultfolder"] + "/codingPotential/mergeCodingPotential.log"
    params:
        outdir = config["resultfolder"] + "/codingPotential",
    shell:
        '''
        python script/mergeCodingPotential.py -c {input.cpc2} -n {input.cnci} -p {input.plek} --pfam {input.pfam} -o {params.outdir} >{log} 2>&1;
        for i in `cat {params.outdir}/noncoding_union.tlist.txt`;do grep $i {input.candidate_gtf};done >{params.outdir}/noncoding_union.gtf 2>>{log};
        cat {params.outdir}/noncoding_union.gtf {input.equal_gtf} >{params.outdir}/transcript_candidate.gtf 2>>{log};
        # cat <(cut -f9 result/cuffcompare/cufcompF.combined.classcode_equal.gtf|cut -d";" -f1|sort -u) <(cut -f9 result/codingPotential/noncoding_union.gtf|cut -d";" -f1|sort -u)|sort|uniq -d|wc -l
        '''