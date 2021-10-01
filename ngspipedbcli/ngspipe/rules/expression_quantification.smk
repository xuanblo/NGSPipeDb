
quantify_method = config['quantify_method'] # htseqcounts or featurecounts
quantify_outdir = join(config["resultsDir"], "quantify")

quantify_ref = config["genomeAnno"] if config['quantify_ref'] == 'genome_annotation' else join(transcript_assembly_outdir, "merged.gtf")

rule reStringtie_by_stringtie:
    message:
        '''
        stringtie was used to assemble transcript
        '''
    input:
        sorted_bam = join(junction_align_outdir, junction_align_method, 'align', "{sample}", "{sample}.sorted.bam"),
        quantify_ref = quantify_ref,
    output:
        gtf = join(quantify_outdir, 'stringtie', 'restringtie', "{sample}", "{sample}.gtf"),
        expr = join(quantify_outdir, 'stringtie', 'restringtie', "{sample}", "{sample}.tab"),
        gene_count = join(quantify_outdir, 'stringtie', 'restringtie', "{sample}", "{sample}_gene_count_piece.csv")
    threads:
        4
    benchmark:
        join(quantify_outdir, 'stringtie', 'restringtie', "{sample}", "benchmark.txt")
    log:
        join(quantify_outdir, 'stringtie', 'restringtie', "{sample}", "run.log")
    shell:
        '''
        stringtie -p {threads} -G {input.quantify_ref} -o {output.gtf} -A {output.expr} -B -e -l {wildcards.sample} {input.sorted_bam} 1>{log} 2>&1;
        echo -e {wildcards.sample}\"\\t\"{output.gtf} >{output.gene_count} 2>>{log};
        '''

rule countMatrix_by_tximport_stringtie:
    message:
        '''
        tximport stringtie was used to calculate fpkm, tpm, counts
        '''
    input:
        ctab = expand(join(quantify_outdir, 'stringtie', 'restringtie', "{sample}", "t_data.ctab"), sample = SAMPLES),
    output:
        gene_count = join(quantify_outdir, 'stringtie', 'quant', "gene.counts.csv"),
    log:
        join(quantify_outdir, 'stringtie', 'quant', "tximport.log.txt")
    params:
        restringtie_dir = join(quantify_outdir, 'stringtie', 'restringtie')
    benchmark:
        join(quantify_outdir, 'stringtie', 'quant', "benchmark.txt")
    shell:
        '''
        Rscript {snake_dir}/scripts/tximport_stringtie_gene_counts.R --stringtie_quant_outdir {params.restringtie_dir} --conditionfile {new_conditionpath} --output {output.gene_count} 1>{log} 2>&1;
        '''

rule countMatrix_by_stringtie:
    message:
        '''
        stringtie was used to calculate fpkm, tpm, counts
        '''
    input:
        gene_count = expand(join(quantify_outdir, 'stringtie', 'restringtie', "{sample}", "{sample}_gene_count_piece.csv"), sample = SAMPLES),
        expr = expand(join(quantify_outdir, 'stringtie', 'restringtie', "{sample}", "{sample}.tab"), sample = SAMPLES)
    output:
        gene_count = join(quantify_outdir, 'stringtie', 'quant', "gene.counts.csv"),
        transcript_count = join(quantify_outdir, 'stringtie', 'quant', "transcript.csv"),
        transcript_fpkm = join(quantify_outdir, 'stringtie', 'quant', "transcript_fpkm_all_samples.tsv"),
        gene_fpkm = join(quantify_outdir, 'stringtie', 'quant', "gene_fpkm_all_samples.tsv")
    log:
        join(quantify_outdir, 'stringtie', 'quant', "prepDE.log")
    benchmark:
        join(quantify_outdir, 'stringtie', 'quant', "benchmark.txt")
    shell:
        '''
        cat {input.gene_count} >{quantify_outdir}/gtflist.txt 2>{log};
        prepDE.py -i {quantify_outdir}/gtflist.txt -t {output.transcript_count} -g {output.gene_count} -l 150 >>{log} 2>&1;
        sed -i 's#|[^,]*,#,#' {output.gene_count};
        # TPM
        python {snake_dir}/scripts/stringtie_expression_matrix.py -m tpm -s {quantify_outdir}/stringtie/restringtie -t {quantify_outdir}/stringtie/quant/transcript_tpm_all_samples.tsv -g {quantify_outdir}/stringtie/quant/gene_tpm_all_samples.tsv >>{log} 2>&1;
        # FPKM
        python {snake_dir}/scripts/stringtie_expression_matrix.py -m fpkm -s {quantify_outdir}/stringtie/restringtie -t {quantify_outdir}/stringtie/quant/transcript_fpkm_all_samples.tsv -g {quantify_outdir}/stringtie/quant/gene_fpkm_all_samples.tsv >>{log} 2>&1;
        # Normalizition
        '''

ruleorder: countMatrix_by_tximport_stringtie > countMatrix_by_stringtie

if config["junction_align_rna_library"] == "fr-firststrand" or config["junction_align_rna_library"] == "--rna-strandness RF":
    featureCounts_library = "-s 2"
elif config["junction_align_rna_library"] == "fr-secondstrand" or config["junction_align_rna_library"] == "--rna-strandness FR":
    featureCounts_library = "-s 1"
elif config["junction_align_rna_library"] == None:
    featureCounts_library = ""

rule quant_counts_by_feature_counts:
    message:
        '''
        featurecounts
        '''
    input:
        quantify_ref = quantify_ref,
        sorted_bam = expand(join(junction_align_outdir, junction_align_method, 'align', "{sample}", "{sample}.sorted.bam"), sample=SAMPLES)
    output:
        gene_counts = join(quantify_outdir, 'featurecounts', 'quant', "gene.counts.csv"),
        transcript_counts = join(quantify_outdir, 'featurecounts', 'quant', "transcript.counts.csv"),
        gene_counts_tmp = join(quantify_outdir, 'featurecounts', 'quant', "gene.counts.tmp"),
        transcript_counts_tmp = join(quantify_outdir, 'featurecounts', 'quant', "transcript.counts.tmp"),
    params:
        featureCounts_library = featureCounts_library
    threads:
        20
    benchmark:
        join(quantify_outdir, 'featurecounts', 'quant', "benchmark.txt")
    log:
        join(quantify_outdir, 'featurecounts', 'quant', "feature_counts.log.txt")
    shell:
        '''
        featureCounts -g gene_id -t exon -p -a {input.quantify_ref} -o {output.gene_counts_tmp} -s 2 -T {threads} {input.sorted_bam} 1>{log} 2>&1;
        featureCounts -g transcript_id -t exon -p -a {input.quantify_ref} -o {output.transcript_counts_tmp} -s 2 -T {threads} {input.sorted_bam} 1>>{log} 2>&1;
        python {snake_dir}/scripts/featurecounts_reformat.py {output.gene_counts_tmp} {output.gene_counts} 1>>{log} 2>&1;
        python {snake_dir}/scripts/featurecounts_reformat.py {output.transcript_counts_tmp} {output.transcript_counts} 1>>{log} 2>&1;
        '''

if config["junction_align_rna_library"] == "fr-firststrand" or config["junction_align_rna_library"] == "--rna-strandness RF":
    htseq_library = "-s reverse"
elif config["junction_align_rna_library"] == "fr-secondstrand" or config["junction_align_rna_library"] == "--rna-strandness FR":
    htseq_library = "-s yes"
elif config["junction_align_rna_library"] == None:
    htseq_library = "-s no"

rule feature_quant_step1_by_htseq_counts:
    message:
        '''
        htseq counts
        '''
    input:
        quantify_ref = quantify_ref,
        sorted_bam = join(junction_align_outdir, junction_align_method, 'align', "{sample}", "{sample}.sorted.bam")
    output:
        htseq_out = join(quantify_outdir, 'htseqcount', 'quant', "{sample}", "{sample}.htseqcount.txt"),
    params:
        htseq_library = htseq_library
    threads:
        5
    benchmark:
        join(quantify_outdir, 'htseqcount', 'quant', "{sample}", "benchmark.txt")
    log:
        join(quantify_outdir, 'htseqcount', 'quant', "{sample}", "feature_quant.log.txt")
    shell:
        '''
        htseq-count {params.htseq_library} -f bam {input.sorted_bam} {input.quantify_ref} >{output.htseq_out} 2>{log};
        '''

rule merge_counts_by_script:
    message:
        '''
        htseq counts merge
        '''
    input:
        htseq_out = expand(join(quantify_outdir, 'htseqcount', 'quant', "{sample}", "{sample}.htseqcount.txt"), sample=SAMPLES)
    output:
        htseqcounts_merge = join(quantify_outdir, 'htseqcount', 'quant', "gene.counts.csv")
    benchmark:
        join(quantify_outdir, 'htseqcount', 'quant', "benchmark.txt")
    log:
        join(quantify_outdir, 'htseqcount', 'quant', "merge_counts.log.txt")
    shell:
        '''
        python {snake_dir}/scripts/merge_htseqcount.py {input.htseq_out} {output.htseqcounts_merge} 1>{log} 2>&1;
        '''
"""
rule htseqcount_report:
    message:
        '''
        htseq counts
        '''
    input:
        counts_stat = expand(config["resultfolder"] + "/htseqcount/{sample}/{sample}.htseqcount.txt", sample=SAMPLES)
    output:
        htseqcount_report = config["resultfolder"] + "/htseqcount_stat_report/reads_counts_stat.csv"
    shell:
        "python script/reads_distrubution.py {input.counts_stat} {output.htseqcount_report}"
        
rule htseqcount:
    message:
        '''
        htseq counts
        '''
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

"""

rule salmon_index:
    '''
    mashmap=`which mashmap`;
    bash {snake_dir}/scripts/generateDecoyTranscriptome.sh -j {threads} -m $mashmap -a {input.transcriptGtf} -g {input.genomeFa} -t {params.outdir}/transcript.fa -o {params.outdir}/getDecoy 1>>{log} 2>&1;
    '''
    input:
        transcriptGtf = quantify_ref,
        genomeFa = config["genomeFasta"]
    output:
        salmon_index = touch(join(quantify_outdir, 'salmon', 'index', 'index.ok'))
    log:
        join(quantify_outdir, 'salmon', 'index', 'index.log.txt')
    params:
        outdir = join(quantify_outdir, 'salmon', 'index')
    threads: 10
    shell:
        '''
        gffread {input.transcriptGtf} -g {input.genomeFa} -w {params.outdir}/transcript.fa 1>{log} 2>&1;
        bash {snake_dir}/scripts/generateDecoyTranscriptome.sh -j {threads} -a {input.transcriptGtf} -g {input.genomeFa} -t {params.outdir}/transcript.fa -o {params.outdir}/getDecoy 1>>{log} 2>&1;
        salmon index -p {threads} --keepDuplicates -t {params.outdir}/getDecoy/gentrome.fa -i {params.outdir}/transcript -k 31 --decoys {params.outdir}/getDecoy/decoys.txt 1>>{log} 2>&1;
        '''

rule salmon_quant:
    input:
        read1 = join(qc_outdir, qc_method, "{sample}", "{sample}.cleanR1.fq.gz"),
        read2 = join(qc_outdir, qc_method,"{sample}", "{sample}.cleanR2.fq.gz"),
        salmon_index = join(quantify_outdir, 'salmon', 'index', 'index.ok'),
        transcriptGtf = quantify_ref,
    output:
        salmon_transcript_quant = join(quantify_outdir, 'salmon', 'quant', '{sample}', 'quant.sf'),
        salmon_gene_quant = join(quantify_outdir, 'salmon', 'quant', '{sample}', 'quant.genes.sf'),
    log: 
        join(quantify_outdir, 'salmon', 'quant', '{sample}', 'quant.log.txt')
    params:
        outdir = join(quantify_outdir, 'salmon', 'quant', '{sample}'),
        index = join(quantify_outdir, 'salmon', 'index', 'transcript')
    threads: 4
    shell:
        '''
        salmon quant -p 4 -l ISR -i {params.index} -1 {input.read1} -2 {input.read2} -o {params.outdir} --validateMappings --seqBias --gcBias --geneMap {input.transcriptGtf} 1>{log} 2>&1;
        '''

rule merge_salmon_gene_tpm:
    input:
        salmon_gene_quant = expand(join(quantify_outdir, 'salmon', 'quant', '{sample}', 'quant.genes.sf'), sample=SAMPLES)
    output:
        merge_matrix = join(quantify_outdir, 'salmon', 'quant', 'gene.tpm.csv')
    log: join(quantify_outdir, 'salmon', 'quant', 'salmon_merge.log.txt')
    shell:
        '''
        python {snake_dir}/scripts/merge_salmon_transcript_tpm.py {quantify_outdir}/salmon/quant {output.merge_matrix} 1>{log} 2>&1;
        '''

rule salmon_counts_by_tximport:
    input:
        salmon_gene_quant = expand(join(quantify_outdir, 'salmon', 'quant', '{sample}', 'quant.sf'), sample=SAMPLES),
        transcriptGtf = quantify_ref,
    output:
        merge_matrix = join(quantify_outdir, 'salmon', 'quant', 'gene.counts.csv')
    params:
        outdir = join(quantify_outdir, 'salmon', 'quant')
    log: join(quantify_outdir, 'salmon', 'quant', 'salmon_counts_by_tximport.log.txt')
    shell:
        '''
        gffread {input.transcriptGtf} -T -o {params.outdir}/transcript.tmp.gtf;
        python {snake_dir}/scripts/tx2gene.py {params.outdir}/transcript.tmp.gtf >{params.outdir}/tx2gene.tsv;
        Rscript {snake_dir}/scripts/tximport_salmon_gene_counts.R --salmon_quant_outdir {params.outdir} --tx2gene {params.outdir}/tx2gene.tsv --conditionfile {new_conditionpath} --output {output.merge_matrix}
        '''

rule quantification:
    message:
        '''
        quantification
        '''
    input:
        gene_count = join(quantify_outdir, quantify_method, 'quant', "gene.counts.csv"),
    output:
        quantification_ok = touch(join(flag_outdir, 'quantification.ok'))