import os
import sys
import pandas as pd

configfile: "config.yaml"

smpList = pd.read_csv(config["sampleListFile"], index_col=0, header=None)
SAMPLES = list(smpList.index)
CURRENTPATH = os.getcwd()
#SAMPLES = ["FF-1"]

# for test
# SAMPLES = config["testSamples"]

###############
# jatropha curcas project
############

rule all:
    input:
        #multiqc = config["resultfolder"] + "/multiqc/multiqc.ok",
        #counts = config["resultfolder"] + "/counts/transcript.csv",
        #sampleQC = config["resultfolder"] + "/sampleQC/sampleQC.ok",
        #diff = config["resultfolder"] + "/diff/diff.ok"
        #cuffcompare = config["resultfolder"] + "/cuffcompare/cufcompF.combined.gtf",
        #trinityout = config["resultfolder"] + "/trinity/trinity_out.ok",
        #lncRNAstat = config["resultfolder"] + "/lncRNA_stat",
        #feelnc = config["resultfolder"] + "/feelnc/lncRNA_classes.txt",
        #fasta = config["resultfolder"] + "/filtermRNA/cufcompF.combined.uix.fa",
        #tlist = config["resultfolder"] + "/codingPotential/noncoding_union.tlist.txt",
        #transcript_fpkm = config["resultfolder"] + "/assembly_final/transcript_fpkm_all_samples.tsv",
        heatmap = config["resultfolder"] + "/sampleQC/sampleQC.ok",
        diff = config["resultfolder"] + "/diff/diff.ok",
        #bam_stat = expand(config["resultfolder"] + "/mapping_stat/{sample}/{sample}.bam_stat.txt", sample=SAMPLES),
        #bams_stat = config["resultfolder"] + "/stat_report/bam_stat_report/reads_mapping_stat.csv"
        #rev_gc_qual = expand(config["resultfolder"] + "/clean_data_stat/{sample}/{sample}_2.seqkit.gc_qual",sample=SAMPLES),
        #clean_data_report = config["resultfolder"] + "/stat_report/clean_data_stat/reads_product.csv",
        #htseqcount =  expand(config["resultfolder"] + "/htseqcount/{sample}/{sample}.htseqcount.txt", sample=SAMPLES)
        #htseqcount_report = config["resultfolder"] + "/htseqcount_stat_report/reads_counts_stat.csv"
        #feature_counts_final = expand(config["resultfolder"] + "/feature_counts_final/{sample}/{sample}.counts_stat.txt", sample=SAMPLES),
        #htseq_counts_final = expand(config["resultfolder"] + "/htseq_counts_final/{sample}/{sample}.htseqcount.txt", sample=SAMPLES),
        #merge_featurecount = config["resultfolder"] + "/feature_counts_final/featurecount.tsv",
        #merge_htseqcount = config["resultfolder"] + "/htseq_counts_final/htseqcount.tsv"
        #anno = config["resultfolder"] + '/genome_stat/summary.txt'
        #dir = config["resultfolder"] + "/eggnog"
        #gene_clusters = config["resultfolder"] + "/timeSerise/result.txt",
        #tissueSpecified_gene = config["resultfolder"] + "/tissueEnrich/ts.tsv",
        #infer = expand(config["resultfolder"] + "/infer_experiment/{sample}/{sample}.infer.txt", sample=SAMPLES)
        #stringtie = expand(config["resultfolder"] + "/assembly/{sample}/transcript.gtf", sample=SAMPLES)
        #cufflinks = expand(config["resultfolder"] + "/cufflinks/{sample}/{sample}.gtf", sample=SAMPLES)
        #counts_stat = expand(config["resultfolder"] + "/feature_counts_stat/{sample}/{sample}.counts_stat.txt", sample=SAMPLES)
        #mapping_star = expand(config["resultfolder"] + "/mapping_star/{sample}/{sample}.sorted.bam.log", sample=SAMPLES)
        #mapping_hisat2 = expand(config["resultfolder"] + "/mapping/{sample}/{sample}.sorted.bam", sample=SAMPLES)
        #genomeIndex = config["resultfolder"] + "/genomeIndex/index.ok"
        #identifyLnc = config["resultfolder"] + "/filter_out/lncRNA.gtf",


rule clean_data_stat:
    input:
        fwd = config["datafolder"] + "/{sample}_1.clean.fq.gz",
        rev = config["datafolder"] + "/{sample}_2.clean.fq.gz"
    output:
        fwd = config["resultfolder"] + "/clean_data_stat/{sample}/{sample}_1.seqkit.summary",
        rev = config["resultfolder"] + "/clean_data_stat/{sample}/{sample}_2.seqkit.summary",
        fwd_gc_qual = config["resultfolder"] + "/clean_data_stat/{sample}/{sample}_1.seqkit.gc_qual",
        rev_gc_qual = config["resultfolder"] + "/clean_data_stat/{sample}/{sample}_2.seqkit.gc_qual",
    threads: 1
    shell:
        '''
        seqkit stat -j {threads} -a {input.fwd} >{output.fwd};
        seqkit stat -j {threads} -a {input.rev} >{output.rev};
        seqkit fx2tab -n -g -q -j {threads} {input.fwd}|cut -f4-|awk '{{gc=gc+$1;qual=qual+$2;}}END{{print gc/NR"\\t"qual/NR}}' >{output.fwd_gc_qual};
        seqkit fx2tab -n -g -q -j {threads} {input.rev}|cut -f4-|awk '{{gc=gc+$1;qual=qual+$2;}}END{{print gc/NR"\\t"qual/NR}}' >{output.rev_gc_qual};
        '''

rule clean_data_report:
    input:
        fwd = expand(config["resultfolder"] + "/clean_data_stat/{sample}/{sample}_1.seqkit.summary", sample=SAMPLES),
        rev = expand(config["resultfolder"] + "/clean_data_stat/{sample}/{sample}_2.seqkit.summary", sample=SAMPLES),
        fwd_gc_qual = expand(config["resultfolder"] + "/clean_data_stat/{sample}/{sample}_1.seqkit.gc_qual", sample=SAMPLES),
        rev_gc_qual = expand(config["resultfolder"] + "/clean_data_stat/{sample}/{sample}_2.seqkit.gc_qual", sample=SAMPLES)
    output:
        clean_data_report = config["resultfolder"] + "/stat_report/clean_data_stat/reads_product.csv",
    shell:
        "python script/reads_product_stat.py {input.fwd} {input.rev} {input.fwd_gc_qual} {input.rev_gc_qual} {output.clean_data_report}"

rule genome_stat:
    input:
        genome = config['genome'],
        gff = config['gff']
    output:
        anno = config["resultfolder"] + '/genome_stat/summary.txt'
    params:
        perlLibPath = os.path.join(CURRENTPATH, config['gaasPerlLib']),
        gaasBin = config['gaas'],
        outdir = config["resultfolder"] + '/genome_stat'
    log:
        fix = config["resultfolder"] + '/genome_stat/fix.log.txt',
        stat = config["resultfolder"] + '/genome_stat/stat.log.txt',
        statistics = config["resultfolder"] + '/genome_stat/statistics.log.txt',
        dist = config["resultfolder"] + '/genome_stat/dist.log.txt',
    conda: "env/gaas.yaml"
    shell:
        '''
            export PERL5LIB=$PERL5LIB:{params.perlLibPath};
            # fix or repair gff file
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

rule genome_stat_report:
    input:
    output:
    shell:
        '''
        '''

rule RemoveRibosomalRNA:
    input:
        fwd = config["datafolder"] + "/{sample}_1.clean.fq.gz",
        rev = config["datafolder"] + "/{sample}_2.clean.fq.gz"
    output:
        fwd = config["resultfolder"] + "/rRNARemoved/{sample}/{sample}.1.non_rRNA.fq",
        rev = config["resultfolder"] + "/rRNARemoved/{sample}/{sample}.2.non_rRNA.fq",
        interleave = config["resultfolder"] + "/rRNARemoved/{sample}/{sample}.reads_interleaved.fq",
        rRNA = config["resultfolder"] + "/rRNARemoved/{sample}/{sample}.reads_rRNA.fq",
        nonrRNA = config["resultfolder"] + "/rRNARemoved/{sample}/{sample}.reads_non_rRNA.fq",
        pe = config["resultfolder"] + "/rRNARemoved/{sample}/{sample}.reads_non_rRNA.fq.pe",
    params:
        mergedFile = config["resultfolder"] + "/rRNARemoved/{sample}/{sample}.reads_interleaved.fq",
        index = config["rRNAIndex"][0] + "," + config["rRNAIndex"][1],
        name = config["resultfolder"] + "/rRNARemoved/{sample}/{sample}."
    log:
        mergeReads = config["resultfolder"] + "/rRNARemoved/{sample}/step1.mergeReads.log",
        sortmerna = config["resultfolder"] + "/rRNARemoved/{sample}/step2.sortmerna.log",
        extractPairedReads = config["resultfolder"] + "/rRNARemoved/{sample}/step3.extractPairedReads.log",
        unmergeReads = config["resultfolder"] + "/rRNARemoved/{sample}/step4.unmergeReads.log"
    conda: "env/rRNA.yaml"
    threads: 5
    shell:
        #"interleave-reads.py {input.fwd} {input.rev} -o {params.mergedFile} 1>{log.mergeReads} 2>&1;"
        "seqtk mergepe {input.fwd} {input.rev} >{params.mergedFile} 2>{log.mergeReads};"
        "sortmerna --ref {params.index} --reads {params.mergedFile} --fastx --aligned {params.name}reads_rRNA --other {params.name}reads_non_rRNA --log -a {threads} -v 1>{log.sortmerna} 2>&1;"
        #"extract-paired-reads.py {params.name}reads_non_rRNA.fq -p {params.name}reads_non_rRNA.fq.pe -s {params.name}reads_non_rRNA.fq.se 1>{log.extractPairedReads} 2>&1;" #有bug
        "seqtk dropse {params.name}reads_non_rRNA.fq >{params.name}reads_non_rRNA.fq.pe;"
        "unmerge-paired-reads.sh {params.name}reads_non_rRNA.fq.pe {params.name}1.non_rRNA.fq {params.name}2.non_rRNA.fq 1>{log.unmergeReads} 2>&1"

rule QC:
    input:
        fwd = config["resultfolder"] + "/rRNARemoved/{sample}/{sample}.1.non_rRNA.fq",
        rev = config["resultfolder"] + "/rRNARemoved/{sample}/{sample}.2.non_rRNA.fq"
    output:
        config["resultfolder"] + "/qc/{sample}/{sample}.1.non_rRNA_val_1_fastqc.zip",
        config["resultfolder"] + "/qc/{sample}/{sample}.2.non_rRNA_val_2_fastqc.zip",
        config["resultfolder"] + "/qc/{sample}/{sample}.1.non_rRNA_val_1.fq",
        config["resultfolder"] + "/qc/{sample}/{sample}.2.non_rRNA_val_2.fq"
    params:
        outfolder = config["resultfolder"] + "/qc/{sample}"
    log: config["resultfolder"] + "/qc/{sample}/qc.log"
    conda: "env/qc.yaml"
    threads: 5
    shell:
        "trim_galore -j {threads} --fastqc --paired {input.fwd} {input.rev} -o {params.outfolder} 1>{log} 2>&1"

rule multiQCStat:
    input:
        expand(config["resultfolder"] + "/qc/{sample}/{sample}.1.non_rRNA_val_1_fastqc.zip", sample=SAMPLES),
        expand(config["resultfolder"] + "/qc/{sample}/{sample}.2.non_rRNA_val_2_fastqc.zip", sample=SAMPLES),
    output: touch(config["resultfolder"] + "/multiqc/multiqc.ok")
    conda: "env/multiqc.yaml"
    params:
        outfolder = config["resultfolder"]
    log: config["resultfolder"] + "/multiqc/multiqc.log"
    shell:
        "multiqc -f -o {params.outfolder}/multiqc -n multiqc {params.outfolder}/qc/ 1>{log} 2>&1"

# http://rseqc.sourceforge.net/

rule genomeIndex:
    input:
        genome = config["genome"],
        gff = config["gff"]
    output:
        genomeIndex = touch(config["resultfolder"] + "/genomeIndex/index.ok")
    conda: "env/mapping.yaml"
    log: config["resultfolder"] + "/genomeIndex/index.log"
    params:
        outfolder = config["resultfolder"] + "/genomeIndex"
    threads: 5
    shell:
        '''
        hisat2_extract_exons.py {input.gff} > {params.outfolder}/genome.exon;
        hisat2_extract_splice_sites.py {input.gff} > {params.outfolder}/genome.ss;
        hisat2-build -p {threads} --ss {params.outfolder}/genome.ss --exon {params.outfolder}/genome.exon {input.genome} {params.outfolder}/genome 1>{log} 2>&1
        '''

# RSeQC check strand
rule mapping_hisat2:
    input:
        fwd = config["resultfolder"] + "/qc/{sample}/{sample}.1.non_rRNA_val_1.fq",
        rev = config["resultfolder"] + "/qc/{sample}/{sample}.2.non_rRNA_val_2.fq",
        indexOk = config["resultfolder"] + "/genomeIndex/index.ok",
    output:
        mapping_hisat2 = config["resultfolder"] + "/mapping/{sample}/{sample}.sorted.bam"
    log:
        config["resultfolder"] + "/mapping/{sample}/{sample}.sorted.bam.log"
    conda: "env/mapping.yaml"
    params:
        index = config["resultfolder"] + "/genomeIndex/JcChina",
        summary = config["resultfolder"] + "/mapping/{sample}/{sample}.mapping.summary",
        outdir = config["resultfolder"] + "/mapping/{sample}"
    threads: 5
    shell:
        '''
        #hisat2 --rna-strandness RF --no-unal --dta -p {threads} --max-intronlen 5000 --summary-file {params.summary} -x {params.index} -1 {input.fwd} -2 {input.rev} --un-conc-gz {params.outdir}/unmapped.fq.gz 1>{log} 2>&1;
        hisat2 --rna-strandness RF --dta -p {threads} --max-intronlen 6000 --summary-file {params.summary} -x {params.index} -1 {input.fwd} -2 {input.rev} 2>{log} |samtools sort -@ {threads} -o {output} 1>>{log} 2>&1;
        '''


rule mapping_star:
    input:
        fwd = config["resultfolder"] + "/qc/{sample}/{sample}.1.non_rRNA_val_1.fq",
        rev = config["resultfolder"] + "/qc/{sample}/{sample}.2.non_rRNA_val_2.fq",
        fasta = config["genome"],
        gff = config["gff"]
    output:
        mapping_star = config["resultfolder"] + "/mapping_star/{sample}/{sample}.Aligned.out.bam"
    log:
        config["resultfolder"] + "/mapping_star/{sample}/{sample}.sorted.bam.log"
    conda: "env/star.yaml"
    params:
        index = config["resultfolder"] + "/mapping_star/JcChinaRfe",
        summary = config["resultfolder"] + "/mapping_star/{sample}/{sample}.mapping.summary",
        outPrefix = config["resultfolder"] + "/mapping_star/{sample}/{sample}."
    threads: 5
    shell:
        '''
        # http://www.bioinfo-scrounger.com/archives/288
        STAR  --runMode genomeGenerate --genomeDir {params.index} --runThreadN {threads} --genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gff} --sjdbOverhang 149;
        STAR --runThreadN {threads} --genomeDir {params.index} --readFilesIn {input.fwd} {input.rev} --outFileNamePrefix {params.outPrefix} --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN {threads} 1>{log} 2>&1;
        '''

'''
rule mapping_tophat:
    input:
        fwd = config["resultfolder"] + "/qc/{sample}/{sample}.1.non_rRNA_val_1.fq",
        rev = config["resultfolder"] + "/qc/{sample}/{sample}.2.non_rRNA_val_2.fq",
        indexOk = config["resultfolder"] + "/genomeIndex/index.ok"
    output:
        config["resultfolder"] + "/mapping/{sample}/{sample}.sorted.bam"
    log:
        config["resultfolder"] + "/mapping/{sample}/{sample}.sorted.bam.log"
    conda: "env/mapping.yaml"
    params:
        index = config["resultfolder"] + "/genomeIndex/JcChina",
        summary = config["resultfolder"] + "/mapping/{sample}/{sample}.mapping.summary",
        outdir = config["resultfolder"] + "/mapping/{sample}"
    threads: 5
    shell:
        "hisat2 --rna-strandness RF --no-unal --dta -p {threads} --max-intronlen 5000 --summary-file {params.summary} -x {params.index} -1 {input.fwd} -2 {input.rev} --un-conc-gz {params.outdir}/unmapped.fq.gz 2>{log} |samtools sort -@ {threads} -o {output} 1>>{log} 2>&1"

'''

rule infer_experiment:
    input:
        bam = config["resultfolder"] + "/mapping/{sample}/{sample}.sorted.bam",
        bed = "../../Database/JcGenome/GCF_000696525.1_JatCur_1.0/cds.bed"
    output:
        infer = config["resultfolder"] + "/infer_experiment/{sample}/{sample}.infer.txt"
    shell:
        '''
        infer_experiment.py -r {input.bed} -i {input.bam} > {output.infer}
        '''

rule bam_stat:
    input:
        bam = config["resultfolder"] + "/mapping/{sample}/{sample}.sorted.bam",
        gff = config["gtf"],
        genome = config["genome"],
    output:
        bam_stat = config["resultfolder"] + "/mapping_stat/{sample}/{sample}.bam_stat.txt"
    log:
        config["resultfolder"] + "/mapping_stat/{sample}/{sample}.log"
    params:
        outdir = config["resultfolder"] + "/mapping_stat/{sample}"
    conda: "env/RSeQC.yaml"
    threads: 5
    shell:
        '''
        # RSeQC
        bam_stat.py -i {input.bam} >{output.bam_stat} 2>{log};
        # qualimap
        qualimap bamqc -bam {input.bam} --java-mem-size=4G --sequencing-protocol strand-specific-reverse -nt {threads} --collect-overlap-pairs -gff {input.gff} -outdir {params.outdir}/qualimap --output-genome-coverage {params.outdir}/qualimap/genome.coverage.txt
        '''

rule bam_stat_report:
    input:
        bams_stat = expand(config["resultfolder"] + "/mapping_stat/{sample}/{sample}.bam_stat.txt", sample=SAMPLES)
    output:
        bam_stat_report = config["resultfolder"] + "/stat_report/bam_stat_report/reads_mapping_stat.csv"
    log: config["resultfolder"] + "/stat_report/bam_stat_report/reads_mapping_stat.log"
    shell:
        "python script/reads_mapping_stat.py {input.bams_stat} {output.bam_stat_report}"

rule htseqcount:
    input:
        gff = config["mRNA"],
        bam = config["resultfolder"] + "/mapping/{sample}/{sample}.sorted.bam"
    output:
        htseqcount = config["resultfolder"] + "/htseqcount/{sample}/{sample}.htseqcount.txt"
    threads: 2
    log:
        config["resultfolder"] + "/htseqcount/{sample}/{sample}.log"
    shell:
        "htseq-count -s reverse -f bam {input.bam} {input.gff} >{output.htseqcount} 2>{log}"

rule feature_counts:
    input:
        gff = config["gtf"],
        bam = config["resultfolder"] + "/mapping/{sample}/{sample}.sorted.bam"
    output:
        counts_stat = config["resultfolder"] + "/feature_counts_stat/{sample}/{sample}.counts_stat.txt"
    threads: 2
    conda: "env/subread.yaml"
    log:
        config["resultfolder"] + "/feature_counts_stat/{sample}/{sample}.log"
    shell:
        "featureCounts -a {input.gff} -o {output.counts_stat} -s 2 -T {threads} {input.bam} 1>{log} 2>&1"


rule htseqcount_report:
    input:
        counts_stat = expand(config["resultfolder"] + "/htseqcount/{sample}/{sample}.htseqcount.txt", sample=SAMPLES)
    output:
        htseqcount_report = config["resultfolder"] + "/htseqcount_stat_report/reads_counts_stat.csv"
    shell:
        "python script/reads_distrubution.py {input.counts_stat} {output.htseqcount_report}"



'''
rule unmappedReads:
    input:
        config["resultfolder"] + "/mapping/{sample}/{sample}.sorted.bam"
    output:
        unmappedbam = temp(config["resultfolder"] + "/unmapped/{sample}/{sample}.unmapped.bam"),
        unmappedfq1 = config["resultfolder"] + "/unmapped/{sample}/{sample}_1.fq",
        unmappedfq2 = config["resultfolder"] + "/unmapped/{sample}/{sample}_2.fq"
    log:
        unmappedbamlog = config["resultfolder"] + "/unmapped/{sample}/{sample}.unmappedbamlog.log",
        bedtools = config["resultfolder"] + "/unmapped/{sample}/{sample}.bedtools.log"
    params:

    conda: "env/mapping.yaml"
    threads: 5
    shell:
        "samtools view -b -f 4 -@ {threads} {input}|samtools sort -@ {threads} -n - 2>{log.unmappedbamlog}|samtools fixmate - {output.unmappedbam} 1>>{log.unmappedbamlog} 2>&1;"
        "bedtools bamtofastq -i {output.unmappedbam} -fq {output.unmappedfq1} -fq2 {output.unmappedfq2} 1>{log.bedtools} 2>&1;"

rule trinity:
    input:
        unmappedfq1 = expand(config["resultfolder"] + "/unmapped/{sample}/{sample}_1.fq", sample=SAMPLES),
        unmappedfq2 = expand(config["resultfolder"] + "/unmapped/{sample}/{sample}_2.fq", sample=SAMPLES),
    output:
        trinitysample = config["resultfolder"] + "/trinity/trinitysample.txt",
        trinityout = touch(config["resultfolder"] + "/trinity/trinity_out.ok")
    params:
        unmappfolder = config["resultfolder"] + "/unmapped/",
        trinityfolder = config["resultfolder"] + "/trinity/trinity_out"
    conda: "env/trinity.yaml"
    threads: 40
    log: trinitysample = config["resultfolder"] + "/trinity/run_trinity.log",
    shell:
        "python script/sample4trinity.py {params.unmappfolder} >{output.trinitysample};"
        "Trinity --seqType fq --max_memory 110G --output {params.trinityfolder} --samples_file {output.trinitysample} --CPU {threads} 1>{log} 2>&1;"

rule unigenePlusGenome:
    input: 
        "unigene.fa",
        "genome.fa",
    output:
        "Ref.fa"
    shell:
        "cat genome.fa unigene.fa >Ref.fa"

rule refIndex:
    input: "Ref.fa"
    output:
        touch(config["resultfolder"] + "/refIndex/index.ok")
    conda: "env/mapping.yaml"
    log: config["resultfolder"] + "/refIndex/index.log"
    params:
        outfolder = config["resultfolder"] + "/genomeIndex"
    shell:
        "hisat2-build {input} {params.outfolder}/ref 1>{log} 2>&1"
'''

rule stringtie:
    input:
        gff = config["gff"],
        bam = config["resultfolder"] + "/mapping/{sample}/{sample}.sorted.bam",
        gtf = config['gtf']
    output:
        stringtie = config["resultfolder"] + "/assembly/{sample}/transcript.gtf"
    conda: "env/stringtie.yaml"
    log: config["resultfolder"] + "/assembly/{sample}/stringtie.log"
    params:
        outfolder = config["resultfolder"] + "/assembly/{sample}"
    threads: 5
    shell:
        '''
        stringtie -p {threads} -G {input.gff} -o {output.stringtie} {input.bam} 1>{log} 2>&1;
        cuffcompare -r {input.gtf} -o {params.outfolder}/cufcompF {output.stringtie};
        less {params.outfolder}/cufcompF.combined.gtf|perl -ne '/transcript_id "(\S+)";.*class_code "(\S)"/;print "$1\\t$2\\n"'|sort -u|cut -f2|sort|uniq -c|perl -ne 's/^\s+//;print "- $_"' >{params.outfolder}/classcode.stat.txt;
        '''

rule stringtieMerge:
    input: 
        gtfs = expand(config["resultfolder"] + "/assembly/{sample}/transcript.gtf", sample = SAMPLES),
        refGtf = config["gff"]
    output: config["resultfolder"] + "/mergedGtf/merged.gtf"
    conda: "env/stringtie.yaml"
    log: config["resultfolder"] + "/mergedGtf/merged.log"
    params:
        outfolder = config["resultfolder"] + "/mergedGtf"
    shell:
        '''
        stringtie --merge -G {input.refGtf} -o {output} {input.gtfs} 1>{log} 2>&1;
        cuffcompare -r {input.refGtf} -o {params.outfolder}/cufcompF {output};
        less {params.outfolder}/cufcompF.combined.gtf|perl -ne '/transcript_id "(\S+)";.*class_code "(\S)"/;print "$1\\t$2\\n"'|sort -u|cut -f2|sort|uniq -c|perl -ne 's/^\s+//;print "- $_"' >{params.outfolder}/classcode.stat.txt;
        '''

rule assembly_cufflinks:
    input:
        gff = config["gff"],
        bam = config["resultfolder"] + "/mapping/{sample}/{sample}.sorted.bam"
    output:
        cufflinks = config["resultfolder"] + "/cufflinks/{sample}/{sample}.gtf"
    log: config["resultfolder"] + "/cufflinks/{sample}/cufflinks.log"
    threads: 5
    shell:
        "cufflinks -o {output.cufflinks} -p {threads} -g {input.gff} --library-type fr-firststrand {input.bam} 1>{log} 2>&1"


rule cuffcompare:
    input:
     mergedGtf = config["resultfolder"] + "/mergedGtf/merged.gtf",
        refGtf = config["gff"]
    output: config["resultfolder"] + "/cuffcompare/cufcompF.combined.gtf"
    params:
        outfolder = directory(config["resultfolder"] + "/cuffcompare"),
        cuffcompare_gtf = config["resultfolder"] + "/cuffcompare/cufcompF.combined.gtf",
    conda: "env/cuffcompare.yaml"
    shell:
        "cuffcompare -G -r {input.refGtf} -o {params.outfolder}/cufcompF {input.mergedGtf}"
        # le result/cuffcompare/cufcompF.combined.gtf|perl -ne '/transcript_id "(\S+)";.*class_code "(\S)"/;print "$1\t$2\n"'|sort -u|cut -f2|sort|uniq -c|perl -ne 's/^\s+//;print "- $_";'

rule feelnc:
    input: 
        refgtf = config["gff"],
        candidaregtf = config["resultfolder"] + "/cuffcompare/cufcompF.combined.gtf",
        genome = config["genome"]
    output:
        classify = config["resultfolder"] + "/feelnc/lncRNA_classes.txt",
    threads: 40
    log:
        filter_log = config["resultfolder"] + "/feelnc/feelnc.log"
    conda: "env/feelnc.yaml"
    params:
        outdir = config["resultfolder"] + "/feelnc"
    shell:
        "python script/repair_gff.py -i {input.refgtf} -o {params.outdir}/repaired.gff;"
        "awk '$3==\"exon\"||$3==\"mRNA\"' {params.outdir}/repaired.gff|grep -P 'gbkey=mRNA' >{params.outdir}/repaired2.gff;"
        "gffread {params.outdir}/repaired2.gff -F -T -o {params.outdir}/repaired3.gtf;"
        "less {params.outdir}/repaired3.gtf|sed 's/transcript_id/;transcript_id/'|awk -F';' '{{print $1$3\";\"$2}}'|sed 's/ gene_id/gene_id/'|sed 's/;/; /' >{params.outdir}/mRNA.gtf;"
        # filter
        "FEELnc_filter.pl -i {input.candidaregtf} --mRNAfile {params.outdir}/mRNA.gtf --monoex=-1 --size=150 -p {threads} -o {log.filter_log} >{params.outdir}/candidate_lncRNA.gtf 2>{log.filter_log}.filter;"
        # coding potential
        "FEELnc_codpot.pl -i {params.outdir}/candidate_lncRNA.gtf -a {params.outdir}/mRNA.gtf -g {input.genome} --mode=shuffle --outdir {params.outdir}/codpot 2>{log.filter_log}.codpot;"
        # class
        "FEELnc_classifier.pl -i {params.outdir}/codpot/candidate_lncRNA.gtf.lncRNA.gtf -a {params.outdir}/mRNA.gtf -l {log.filter_log}.classifier > {params.outdir}/lncRNA_classes.txt ;"
        "gffread {params.outdir}/codpot/candidate_lncRNA.gtf.lncRNA.gtf -g {input.genome} -w {params.outdir}/candidate_lncRNA.fa;"
        #"FEELnc_pipeline.sh --candidate=result/cuffcompare/cufcompF.combined.gtf --reference=result/feelnc/repaired.gtf --genome=../../Database/JcGenome/GCF_000696525.1_JatCur_1.0/GCF_000696525.1_JatCur_1.0_genomic.fna --outname=FEELnc_pipeline --outdir=result/feelnc >{log} 2>&1"

# filter lncRNA by class_code (uix)
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
'''

rule cpc2:
    input: config["resultfolder"] + "/feelnc/candidate_lncRNA.fa"
    output: config["resultfolder"] + "/codingPotential/cpc2.out.txt"
    conda: "env/cnci.yaml"
    log: config["resultfolder"] + "/codingPotential/cpc2.log"
    shell:
        "python ../../Software/CPC2-beta/bin/CPC2.py -i {input} -o {output} >{log} 2>&1"

rule cnci:
    input: config["resultfolder"] + "/feelnc/candidate_lncRNA.fa"
    output: config["resultfolder"] + "/codingPotential/CNCI.index"
    conda: "env/cnci.yaml"
    params:
        outdir = config["resultfolder"] + "/codingPotential"
    threads: 40
    log: config["resultfolder"] + "/codingPotential/cnci.log"
    shell:
        "cd ../../Software;"
        "python CNCI_package/CNCI.py -f ../Analysis/RunAll/{input} -o ../Analysis/RunAll/{params.outdir} -m pl -p {threads} >../Analysis/RunAll/{log} 2>&1"

# plek change sys.exit(1) to sys.exit(0)
rule plek:
    input: config["resultfolder"] + "/feelnc/candidate_lncRNA.fa"
    output: config["resultfolder"] + "/codingPotential/PLEK_out"
    threads: 40
    log: config["resultfolder"] + "/codingPotential/PLEK.log"
    conda: "env/plek.yaml"
    shell:
        "python ../../Software/PLEK.1.2/PLEK.py -fasta {input} -out {output} -thread {threads} >{log} 2>&1" # line 407 sys.exit(0)

# https://pypi.org/project/gff3tool/

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
        "python script/mergeCodingPotential.py -c {input.cpc2} -n {input.cnci} -p {input.plek} --pfam {input.pfam} -o {params.outdir} >{log} 2>&1;"
        "for i in `cat {params.outdir}/noncoding_union.tlist.txt <(cut -f3 {params.feelncdir}/lncRNA_classes.txt|grep -v 'lncRNA_transcript'|sort -u)|sort|uniq -d`;do grep $i {params.feelncdir}/candidate_lncRNA.gtf;done >{params.outdir}/noncoding_union.gtf;"
        "cat {params.outdir}/noncoding_union.gtf {input.mRNA_gtf} >{params.outdir}/transcript_candidate.gtf;"

rule reStringtie:
    input:
        transcript_candidate_gtf = config["resultfolder"] + "/codingPotential/transcript_candidate.gtf",
        bam = config["resultfolder"] + "/mapping/{sample}/{sample}.sorted.bam"
    output:
        gtf = config["resultfolder"] + "/quant/{sample}/{sample}.gtf",
        expr = config["resultfolder"] + "/quant/{sample}/{sample}.tab",
        gene_count = config["resultfolder"] + "/quant/{sample}/{sample}_gene_count_piece.csv",
    conda: "env/stringtie.yaml"
    threads: 4
    params:
        outdir = config["resultfolder"] + "/quant"
    shell:
        "stringtie -p {threads} -G {input.transcript_candidate_gtf} -o {output.gtf} -A {output.expr} -B -e -l {wildcards.sample} {input.bam};"
        "echo -e {wildcards.sample}\"\\t\"{output.gtf} >{output.gene_count}"

rule countMatrix:
    input:
        gene_count = expand(config["resultfolder"] + "/quant/{sample}/{sample}_gene_count_piece.csv", sample = SAMPLES),
        expr = expand(config["resultfolder"] + "/quant/{sample}/{sample}.tab", sample = SAMPLES)
    output:
        gene_count = config["resultfolder"] + "/counts/gene.csv",
        transcript_count = config["resultfolder"] + "/counts/transcript.csv",
        transcript_fpkm = config["resultfolder"] + "/counts/transcript_fpkm_all_samples.tsv",
        gene_fpkm = config["resultfolder"] + "/counts/gene_fpkm_all_samples.tsv"
    params:
        outfolder = config["resultfolder"] + "/counts",
        stringtiefolder = config["resultfolder"] + "/quant"
    log: config["resultfolder"] + "/counts/prepDE.log"
    shell:
        "cat {input.gene_count} >{params.outfolder}/gtflist.txt;"
        "prepDE.py -i {params.outfolder}/gtflist.txt -t {output.transcript_count} -g {output.gene_count} -l 150 >{log} 2>&1;"
        # TPM
        "python script/stringtie_expression_matrix.py -m tpm -s {params.stringtiefolder} -t {params.outfolder}/transcript_tpm_all_samples.tsv -g {params.outfolder}/gene_tpm_all_samples.tsv;"
        # FPKM
        "python script/stringtie_expression_matrix.py -m fpkm -s {params.stringtiefolder} -t {params.outfolder}/transcript_fpkm_all_samples.tsv -g {params.outfolder}/gene_fpkm_all_samples.tsv;"
        # Normalizition

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

rule stat_lncRNA:
    input:
        gtf = config["resultfolder"] + "/filter_out/lncRNA.gtf",
    output:
        directory(config["resultfolder"] + "/lncRNA_stat")
    conda: "env/statlncRNA.yaml"
    log: config["resultfolder"] + "/lncRNA_stat/lncRNA_stat.log"
    shell:
        "python script/gff_stat_visualizition_v1.py -i {input.gtf} -o {output} >{log} 2>&1;"
        "python script/filter_lncRNA_stat.py >xxx"

rule final_data:
    input:
        mRNA_gtf = config["resultfolder"] + "/feelnc/mRNA.gtf",
        lncRNA_gtf = config["resultfolder"] + "/filter_out/lncRNA.gtf",
        gene_count = config["resultfolder"] + "/counts/gene.csv",
        transcript_count = config["resultfolder"] + "/counts/transcript.csv",
        transcript_fpkm = config["resultfolder"] + "/counts/transcript_fpkm_all_samples.tsv",
        gene_fpkm = config["resultfolder"] + "/counts/gene_fpkm_all_samples.tsv",
        genome = config["genome"]
    output:
        mRNA_fa = config["resultfolder"] + "/assembly_final/mRNA.fa",
        lncRNA_fa = config["resultfolder"] + "/assembly_final/lncRNA.fa",
        gene_count = config["resultfolder"] + "/assembly_final/gene.tsv",
        transcript_count = config["resultfolder"] + "/assembly_final/transcript.tsv",
        transcript_fpkm = config["resultfolder"] + "/assembly_final/transcript_fpkm_all_samples.tsv",
        gene_fpkm = config["resultfolder"] + "/assembly_final/gene_fpkm_all_samples.tsv",
        mRNA_gtf = config["resultfolder"] + "/assembly_final/mRNA.gtf",
        lncRNA_gtf = config["resultfolder"] + "/assembly_final/lncRNA.gtf",
        transcript_gtf = config["resultfolder"] + "/assembly_final/transcript.gtf",
        gene_count_norm = config["resultfolder"] + "/assembly_final/gene.counts.norm.tsv",
        gene_fpkm_norm = config["resultfolder"] + "/assembly_final/gene.fpkm.norm.tsv",
        transcript_count_norm = config["resultfolder"] + "/assembly_final/transcript.counts.norm.tsv",
        transcript_fpkm_norm = config["resultfolder"] + "/assembly_final/transcript.fpkm.norm.tsv",
    log:
        outdir = config["resultfolder"] + "/assembly_final"
    shell:
        '''
        cp {input.mRNA_gtf} {output.mRNA_gtf};
        cp {input.lncRNA_gtf} {output.lncRNA_gtf};
        cat {input.mRNA_gtf} {output.lncRNA_gtf} >{output.transcript_gtf};
        gffread {output.mRNA_gtf} -g {input.genome} -w {output.mRNA_fa};
        gffread {output.lncRNA_gtf} -g {input.genome} -w -T -o {output.lncRNA_fa};
        python script/extract_exp.py -i {output.transcript_gtf} -e {input.gene_count} -o {output.gene_count};
        python script/extract_exp.py -i {output.transcript_gtf} -e {input.gene_fpkm} -o {output.gene_fpkm};
        python script/extract_exp.py -i {output.transcript_gtf} -e {input.transcript_count} -o {output.transcript_count};
        python script/extract_exp.py -i {output.transcript_gtf} -e {input.transcript_fpkm} -o {output.transcript_fpkm};
        perl ../../Software/RNAseq2Coexpnetwork/Final/data_normalization.pl --gene_matrix={output.gene_count} --output_file={output.gene_count_norm} --normalization=deseq --outputstat >{log.outdir}/gene_count_norm.log 2>&1;
        perl ../../Software/RNAseq2Coexpnetwork/Final/data_normalization.pl --gene_matrix={output.gene_fpkm} --output_file={output.gene_fpkm_norm} --normalization=deseq --outputstat >{log.outdir}/gene_fpkm_norm.log 2>&1;
        perl ../../Software/RNAseq2Coexpnetwork/Final/data_normalization.pl --gene_matrix={output.transcript_count} --output_file={output.transcript_count_norm} --normalization=deseq --outputstat >{log.outdir}/transcript_count_norm.log 2>&1;
        perl ../../Software/RNAseq2Coexpnetwork/Final/data_normalization.pl --gene_matrix={output.transcript_fpkm} --output_file={output.transcript_fpkm_norm} --normalization=deseq --outputstat >{log.outdir}/transcript_fpkm_norm.log 2>&1;
        '''

rule htseq_counts_final:
    input:
        transcript_gtf = config["resultfolder"] + "/assembly_final/transcript.gtf",
        bam = config["resultfolder"] + "/mapping/{sample}/{sample}.sorted.bam"
    output:
        htseq_counts_final = config["resultfolder"] + "/htseq_counts_final/{sample}/{sample}.htseqcount.txt"
    threads: 2
    log:
        config["resultfolder"] + "/htseq_counts_final/{sample}/{sample}.log"
    shell:
        "htseq-count -s reverse -f bam {input.bam} {input.transcript_gtf} >{output.htseq_counts_final} 2>{log}"

rule feature_counts_final:
    input:
        transcript_gtf = config["resultfolder"] + "/assembly_final/transcript.gtf",
        bam = config["resultfolder"] + "/mapping/{sample}/{sample}.sorted.bam"
    output:
        feature_counts_final = config["resultfolder"] + "/feature_counts_final/{sample}/{sample}.counts_stat.txt"
    threads: 2
    conda: "env/subread.yaml"
    log:
        config["resultfolder"] + "/feature_counts_final/{sample}/{sample}.log"
    shell:
        "featureCounts -a {input.transcript_gtf} -p -o {output.feature_counts_final} -s 2 -T {threads} {input.bam} 1>{log} 2>&1"

rule merge_featurecount:
    input: expand(config["resultfolder"] + "/feature_counts_final/{sample}/{sample}.counts_stat.txt", sample=SAMPLES)
    output:
        merge_featurecount = config["resultfolder"] + "/feature_counts_final/featurecount.tsv"
    shell:
        '''
        python script/merge_featureCount.py {input} {output}
        '''

rule merge_htseqcount:
    input: expand(config["resultfolder"] + "/htseq_counts_final/{sample}/{sample}.htseqcount.txt", sample=SAMPLES)
    output:
        merge_htseqcount = config["resultfolder"] + "/htseq_counts_final/htseqcount.tsv"
    shell:
        '''
        python script/merge_htseqcount.py {input} {output}
        '''

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

rule QCsamples:
    input:
        gene_count = config["resultfolder"] + "/assembly_final/gene.counts.norm.tsv",
    output:
        ok = touch(config["resultfolder"] + "/sampleQC/sampleQC.ok")
    params:
        countFolder = config["resultfolder"] + "/sampleQC"
    conda: "env/trinity.yaml"
    log: config["resultfolder"] + "/sampleQC/QCsamples.log"
    shell:
        "cd {params.countFolder};"
        # generate a barplot that sums frag counts per replicate across all samples.
        "PtR --matrix ../../result/assembly_final/gene.tsv --samples ../../samples.txt --log2 --barplot_sum_counts >../../{log} 2>&1;"
        # generate a sample correlation matrix plot
        "PtR --matrix ../../result/assembly_final/gene.tsv --samples ../../samples.txt --log2 --CPM --sample_cor_matrix >../../{log}1 2>&1;"
        # PCA
        "PtR --matrix ../../result/assembly_final/gene.tsv --samples ../../samples.txt --log2 --CPM --prin_comp 2 >../../{log}2 2>&1"

rule DEanalysis:
    input:
        gene_count = config["resultfolder"] + "/assembly_final/gene.tsv",
    output:
        ok = touch(config["resultfolder"] + "/diff/diff.ok")
    conda: "env/trinity.yaml"
    log: config["resultfolder"] + "/diff/diff.log"
    params:
        diffFolder = config["resultfolder"] + "/diff"
    shell:
        "run_DE_analysis.pl --matrix {input.gene_count} --method DESeq2 --samples_file samples.txt --output {params.diffFolder} >{log} 2>&1;"
        # Extracting and clustering differentially expressed transcripts
        "cd {params.diffFolder};"
        "analyze_diff_expr.pl --matrix ../../result/assembly_final/gene.tsv --samples ../../samples.txt >../../{log}1 2>&1;"
        # Automatically Partitioning Genes into Expression Clusters
        "define_clusters_by_cutting_tree.pl -R diffExpr.P0.001_C2.matrix.RData --Ptree 60 >../../{log}2 2>&1;"


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

rule coexp:
    input: 
        gene_count = config["resultfolder"] + "/assembly_final/gene.tsv",
    output:
        normalized_count_matrix = config["resultfolder"] + "/coexp/gene_counts_matrix.tsv",
        coexp_network = config["resultfolder"] + "/coexp/co-exp.network",
        stat = config["resultfolder"] + "/coexp/stat.txt",
    log:
        normalized_log = config["resultfolder"] + "/coexp/normalized.log",
        network_construction_log = config["resultfolder"] + "/coexp/network_construction.log",
        network_stat_log = config["resultfolder"] + "/coexp/network_stat.log",
    threads: 40
    shell:
        "perl network_construction.pl --gene_matrix={output.normalized_count_matrix} --output_network={output.coexp_network} --log2x --filter_mean=FLOAT(+) --filter_sd=FLOAT(+) --cc_method=pcc --cutoff_adjpv=0.05 --cutoff_ccv=0 --cpucore={threads} --signed --outputstat >{log.normalized_log} 2>&1;"
        "perl network_stat.pl --network={output.coexp_network} --degree_stat={output.stat} >{log.normalized_log} 2>&1;"

'''
rule tar:
    input:
    output: "fungi_transcriptome_analysis.tar.gz"
    shell:
        "tar -zcvf result/fungi_transcriptome_analysis.tar.gz result/sampleQC result/diff result/counts result/eggnog result/transdecoder result/gtf2statdir/"
'''
