rule gffcompare:
    '''
    deal with novel isoform
    '''
    input:
        mergedGtf = join(transcript_assembly_outdir, "merged.gtf"),
        refGtf = config["genomeAnno"],
        genome = config["genomeFasta"]
    output:
        compared_gtf = join(identify_lncRNA_outdir, 'gffcompare', 'gffcompF.annotated.gtf'),
        remove_class_equal_fa = join(identify_lncRNA_outdir, 'gffcompare', 'gffcompF.annotated.classcode_except_equal.fa'),
        except_queal_gtf = join(identify_lncRNA_outdir, 'gffcompare', 'gffcompF.annotated.classcode_except_equal.gtf'),
    params:
        outfolder = join(identify_lncRNA_outdir, 'gffcompare')
    log: join(identify_lncRNA_outdir, 'gffcompare', 'run.log.txt')
    shell:
        '''
        gffcompare -r {input.refGtf} -o {params.outfolder}/gffcompF {input.mergedGtf} 1>{log} 2>&1;
        # remove no strand transcript id list
        less {output.compared_gtf}|awk '$7=="."'|cut -f9|cut -d" " -f2|sed 's/[";]//g'|sort -u >{params.outfolder}/gffcompF.annotated.noStrand.tlist 2>>{log};
        # class code not equal to "=" transcript id
        less {output.compared_gtf}|awk '$3=="transcript"'|grep -v -P 'class_code "="'|perl -ne '/transcript_id "(\S+)";/;print "$1\n"'|sort -u >{params.outfolder}/gffcompF.annotated.classcode_except_equal.tlist 2>>{log};
        # gtf remove class code "="
        perl {snake_dir}/scripts/fish_gtf_by_tid.pl {params.outfolder}/gffcompF.annotated.classcode_except_equal.tlist {output.compared_gtf} >{params.outfolder}/gffcompF.annotated.classcode_except_equal.gtf 2>>{log};
        # class code equal to "=" transcript id
        less {output.compared_gtf}|awk '$3=="transcript"'|grep -P 'class_code "="'|perl -ne '/transcript_id "(\S+)";/;print "$1\n"'|sort -u >{params.outfolder}/gffcompF.annotated.classcode_equal.tlist 2>>{log};
        # gtf class code "="
        perl {snake_dir}/scripts/fish_gtf_by_tid.pl {params.outfolder}/gffcompF.annotated.classcode_equal.tlist {output.compared_gtf} >{params.outfolder}/gffcompF.annotated.classcode_equal.gtf 2>>{log};
        # novol transcript fasta
        gffread {params.outfolder}/gffcompF.annotated.classcode_except_equal.gtf -g {input.genome} -w {params.outfolder}/gffcompF.annotated.classcode_except_equal.fa 1>>{log} 2>&1;
        # all transcript fasta
        gffread {output.compared_gtf} -g {input.genome} -w {params.outfolder}/gffcompF.annotated.fa  1>>{log} 2>&1;
        # statistic class code
        less {output.compared_gtf}|perl -ne '/transcript_id "(\S+)";.*class_code "(\S)"/;print "$1\\t$2\\n"'|sort -u|cut -f2|sort|uniq -c|perl -ne 's/^\s+//;print "- $_"' >{params.outfolder}/classcode.stat.txt 2>>{log};
        '''

rule feelnc:
    input:
        except_queal_gtf = join(identify_lncRNA_outdir, 'gffcompare', 'gffcompF.annotated.classcode_except_equal.gtf'),
        mRNAgtf = config["genomeAnno"],
        genomeFa = config["genomeFasta"],
    output:
        feelnc_lncRNA_gtf = join(identify_lncRNA_outdir, 'feelnc', 'codpot', 'ngspipe.codpot.lncRNA.gtf'),
        feelnc_mRNA_gtf = join(identify_lncRNA_outdir, 'feelnc', 'codpot', 'ngspipe.codpot.mRNA.gtf'),
    params:
        outdir = join(identify_lncRNA_outdir, 'feelnc'),
        prefix = 'ngspipe',
    threads: 40
    log:
        filter = join(identify_lncRNA_outdir, 'feelnc', 'run.log.txt'),
    shell:
        '''
        export FEELNCPATH=~/miniconda3/pkgs/feelnc-0.2-pl526_0;
        # 1. Filter
        # FEELnc_filter.pl -i transcript_chr38.gtf -a annotation_chr38.gtf -b transcript_biotype=protein_coding > candidate_lncRNA.gtf
        # 2. Coding_Potential
        #FEELnc_codpot.pl -i {input.except_queal_gtf} -a {input.mRNAgtf} -g {input.genomeFa} --outdir={params.outdir} -o {params.prefix} --mode intergenic 1>{log} 2>&1;
        # 3. Classifier
        #FEELnc_classifier.pl -i feelnc_codpot_out/candidate_lncRNA.gtf.lncRNA.gtf -a annotation_chr38.gtf > candidate_lncRNA_classes.txt
        # or use pipeline
        FEELnc_pipeline.sh --candidate={input.except_queal_gtf} --reference={input.mRNAgtf} --genome={input.genomeFa} --outname={params.prefix} --outdir={params.outdir} 1>{log} 2>&1;
        '''

rule filter_lncRNA:
    input:
        feelnc_lncRNA_gtf = join(identify_lncRNA_outdir, 'feelnc', 'codpot', 'ngspipe.codpot.lncRNA.gtf'),
    output:
        filted_lncRNA_gtf = join(identify_lncRNA_outdir, 'filted_lncRNA', 'filted_lncRNA.gtf'),
    params:
        filter_out = join(identify_lncRNA_outdir, 'filted_lncRNA'),
    log: join(identify_lncRNA_outdir, 'filted_lncRNA', 'run.log.txt'),
    shell:
        '''
        python {snake_dir}/scripts/filter_exon_step1.py -i {input.feelnc_lncRNA_gtf} -n 1 -o {params.filter_out}/step1_exon_filter.gtf 1>{log} 2>&1;
        python {snake_dir}/scripts/filter_transcript_length_step2.py -i {params.filter_out}/step1_exon_filter.gtf -l 200 -o {params.filter_out}/step2_length_filter.gtf 1>>{log} 2>&1;
        ln -s {params.filter_out}/step2_length_filter.gtf {output.filted_lncRNA_gtf} 1>>{log} 2>&1;
        '''

rule merge_ref_and_lncRNA:
    input:
        filted_lncRNA_gtf = join(identify_lncRNA_outdir, 'filted_lncRNA', 'filted_lncRNA.gtf'),
        mRNAgtf = config["genomeAnno"],
    output:
        final_gtf = join(identify_lncRNA_outdir, 'final', 'final.gtf'),
    log: join(identify_lncRNA_outdir, 'final', 'run.log.txt')
    shell:
        '''
        cat {input.filted_lncRNA_gtf} {input.mRNAgtf}|bedtools sort -i - >{output.final_gtf} 2>{log};
        '''

config['genomeAnno'] = join(identify_lncRNA_outdir, 'final', 'final.gtf')

"""

# 从有参考基因组的转录结果GTF文件预测编码区域
rule transdecoder:
    input:
        lncRNAtlist = config["resultfolder"] + "/feelnc/lncRNA_final.tlist",
        except_queal_gtf = config["resultfolder"] + "/gffcompare/gffcompF.annotated.classcode_except_equal.gtf",
        genomePep = config['genome_pep'],
        wgj_pep = config['wgj_pep'],
        genomeFa = config['genome']
    output:
        proteinFa = config["resultfolder"] + "/transdecoder/protein.fa",
        TUCP_list = config["resultfolder"] + "/transdecoder/TUCP_list.tsv",
    conda: 'env/transdecoder.yaml'
    params:
        LongOrfs_dir = config["resultfolder"] + "/transdecoder/LongOrfs_dir",
        outdir = config["resultfolder"] + "/transdecoder",
    log: "run.log",
    shell:
        '''
        cat {input.lncRNAtlist} <(cut -f9 {input.except_queal_gtf}|cut -d";" -f1|cut -d" " -f2|sed 's/"//g'|sort -u)|sort|uniq -u >{params.LongOrfs_dir}/candidate_new_mRNA.tlist;
        perl script/fish_gtf_by_tid.pl {params.LongOrfs_dir}/candidate_new_mRNA.tlist {input.except_queal_gtf} >{params.LongOrfs_dir}/candidate_new_mRNA.gtf;
        gffread {params.LongOrfs_dir}/candidate_new_mRNA.gtf -g {input.genomeFa} -w {params.LongOrfs_dir}/protein_candidates.fa;
        TransDecoder.LongOrfs -S -t {params.LongOrfs_dir}/protein_candidates.fa --output_dir {params.LongOrfs_dir} >{log} 2>&1;
        makeblastdb -in {input.wgj_pep} -dbtype prot -out {params.outdir}/wgj_pep/index;
        blastp -query {params.LongOrfs_dir}/longest_orfs.pep -db {params.outdir}/wgj_pep/index -max_target_seqs 1 -outfmt 6 -evalue 1e-10 -num_threads 40 -out {params.outdir}/blastp/blastp.txt;
        #hmmscan --cpu 8 --domtblout pfam.domtblout /path/to/Pfam-A.hmm transdecoder_dir/longest_orfs.pep;
        cd result/transdecoder/predict;
        TransDecoder.Predict -t ../../../{params.LongOrfs_dir}/protein_candidates.fa --output_dir ../../../{params.LongOrfs_dir} --retain_blastp_hits ../../../{params.outdir}/blastp/blastp.txt >>{log} 2>&1;
        cd ../../../;
        less result/transdecoder/predict/protein_candidates.fa.transdecoder.bed |grep -v "track"|cut -f1|sort -u >result/transdecoder/new_mRNA.tlist;
        perl script/fish_gtf_by_tid.pl result/transdecoder/new_mRNA.tlist {input.except_queal_gtf} >{params.outdir}/mRNA_final.gtf;
        cat result/transdecoder/new_mRNA.tlist result/transdecoder/LongOrfs_dir/candidate_new_mRNA.tlist |sort|uniq -u >result/transdecoder/tucp.tlist;
        perl script/fish_gtf_by_tid.pl result/transdecoder/tucp.tlist {input.except_queal_gtf} >{params.outdir}/tucp_final.gtf;
        cat <(cut -d"-" -f1 result/gffcompare/gffcompF.annotated.classcode_equal.tlist) <(grep ">" ../../Database/JcGenome/cms/newgenomejc.proteins.fasta|cut -d" " -f1|sed 's/>//'|sort -u)|sort|uniq -d >result/transdecoder/gffcompF.annotated.classcode_equal.mRNA.tlist;
        perl script/fishInWinter.pl -bf table -ff fasta -gene result/transdecoder/gffcompF.annotated.classcode_equal.mRNA.tlist ../../Database/JcGenome/cms/newgenomejc.proteins.fasta >result/transdecoder/gffcompF.annotated.classcode_equal.mRNA.fa;
        cat result/transdecoder/gffcompF.annotated.classcode_equal.mRNA.fa result/transdecoder/predict/protein_candidates.fa.transdecoder.pep|sed 's/*$//' >result/transdecoder/all_protein.fa;
        '''

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
"""