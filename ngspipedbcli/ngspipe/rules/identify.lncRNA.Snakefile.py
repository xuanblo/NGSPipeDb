rule filter_exon_exp:
    input:
        feelnc = config["resultfolder"] + "/feelnc/lncRNA_classes.txt",
        merged_gtf = config["resultfolder"] + "/cuffcompare/cufcompF.combined.gtf",
        exp = config["resultfolder"] + "/counts/transcript_fpkm_all_samples.tsv",
    output:
        filtered_lncRNA = config["resultfolder"] + "/filter_exon_exp/lncRNA.gtf"
    log:
        config["resultfolder"] + "/filter_exon_exp/filter_exon_exp.log"
    shell:
        '''
        python script/filter_exon_exp.py {input.feelnc} {input.merged_gtf} {input.exp} >{output.filtered_lncRNA} 2>{log};
        '''

rule identifyLnc:
    input:
        gtf = config["resultfolder"] + "/codingPotential/noncoding_union.gtf",
        tlist = config["resultfolder"] + "/codingPotential/noncoding_union.tlist.txt",
        fpkm_matrix = config["resultfolder"] + "/counts/transcript_fpkm_all_samples.tsv",
        refGtf = config["gtf"]
    output:
        identifyLnc = config["resultfolder"] + "/filter_out/lncRNA.gtf",
    params:
        filter_out = config["resultfolder"] + "/filter_out",
    conda: "env/matplotlib.yaml"
    shell:
        '''
        mkdir -p {params.filter_out};
        # step1
        python script/filter_exon_step1.py -i {input.gtf} -n 1 -o {params.filter_out}/step1_exon_filter.gtf;
        cuffcompare -r {input.refGtf} -o {params.filter_out}/step1 {params.filter_out}/step1_exon_filter.gtf;
        echo "gene:" `less {params.filter_out}/step1.combined.gtf |cut -f9|cut -d";" -f1|sort -u|wc -l` "transcript:" `less {params.filter_out}/step1.combined.gtf |cut -f9|cut -d";" -f2|sort -u|wc -l` >{params.filter_out}/step1.featureNum.stat.txt
        less {params.filter_out}/step1.combined.gtf|perl -ne '/transcript_id "(\S+)";.*class_code "(\S)"/;print "$1\\t$2\\n"'|sort -u|cut -f2|sort|uniq -c|perl -ne 's/^\s+//;print "- $_"' >{params.filter_out}/step1.classcode.stat.txt;
        # step2
        python script/filter_transcript_length_step2.py -i {params.filter_out}/step1_exon_filter.gtf -l 200 -o {params.filter_out}/step2_length_filter.gtf;
        cuffcompare -r {input.refGtf} -o {params.filter_out}/step2 {params.filter_out}/step2_length_filter.gtf;
        echo "gene:" `less {params.filter_out}/step2.combined.gtf |cut -f9|cut -d";" -f1|sort -u|wc -l` "transcript:" `less {params.filter_out}/step2.combined.gtf |cut -f9|cut -d";" -f2|sort -u|wc -l` >{params.filter_out}/step2.featureNum.stat.txt
        less {params.filter_out}/step2.combined.gtf|perl -ne '/transcript_id "(\S+)";.*class_code "(\S)"/;print "$1\\t$2\\n"'|sort -u|cut -f2|sort|uniq -c|perl -ne 's/^\s+//;print "- $_"' >{params.filter_out}/step2.classcode.stat.txt;
        # step3
        python script/filter_low_expression_transcript_step3.py -i {params.filter_out}/step2_length_filter.gtf -m {input.fpkm_matrix} -f 0.5 -o {params.filter_out}/step3_exp_filter.gtf;
        cuffcompare -r {input.refGtf} -o {params.filter_out}/step3 {params.filter_out}/step3_exp_filter.gtf;
        echo "gene:" `less {params.filter_out}/step3.combined.gtf |cut -f9|cut -d";" -f1|sort -u|wc -l` "transcript:" `less {params.filter_out}/step3.combined.gtf |cut -f9|cut -d";" -f2|sort -u|wc -l` >{params.filter_out}/step3.featureNum.stat.txt
        less {params.filter_out}/step3.combined.gtf|perl -ne '/transcript_id "(\S+)";.*class_code "(\S)"/;print "$1\\t$2\\n"'|sort -u|cut -f2|sort|uniq -c|perl -ne 's/^\s+//;print "- $_"' >{params.filter_out}/step3.classcode.stat.txt;
        # step4
        python script/filter_coding_potential_step4.py -i {params.filter_out}/step3_exp_filter.gtf -t {input.tlist} -o {params.filter_out}/lncRNA.gtf;
        cuffcompare -r {input.refGtf} -o {params.filter_out}/step4 {params.filter_out}/lncRNA.gtf;
        echo "gene:" `less {params.filter_out}/step4.combined.gtf |cut -f9|cut -d";" -f1|sort -u|wc -l` "transcript:" `less {params.filter_out}/step4.combined.gtf |cut -f9|cut -d";" -f2|sort -u|wc -l` >{params.filter_out}/step4.featureNum.stat.txt
        less {params.filter_out}/step4.combined.gtf|perl -ne '/transcript_id "(\S+)";.*class_code "(\S)"/;print "$1\\t$2\\n"'|sort -u|cut -f2|sort|uniq -c|perl -ne 's/^\s+//;print "- $_"' >{params.filter_out}/step4.classcode.stat.txt;
        '''
