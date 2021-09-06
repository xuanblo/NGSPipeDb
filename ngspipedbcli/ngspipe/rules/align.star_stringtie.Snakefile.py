# ----------------------------------------------------------------------------
# step1
step1_outdir = join(config["result_dir"], "align.hisat2_stringtie", \
    "step1_index_genome_by_hisat2_build")
step1_dbname = "genome"
# ----------------------------------------------------------------------------
rule align_hisat2_stringtie_step1_index_genome_by_hisat2_build:
    message:
        '''
        ------------------------------
        hisat2-build
        ------------------------------
        '''
    input:
        genome_fasta = config["genome_fasta"],
        genome_gtf = config["genome_gtf"],
    output:
        touch(join(step1_outdir, "index.ok"))
    params:
        dbPrefix = join(step1_outdir, step1_dbname),
        outdir = step1_outdir
    conda:
        py3env
    benchmark:
        join(step1_outdir, "benchmark.txt")
    log:
        join(step1_outdir, "index.log")
    threads:
        5
    shell:
        '''
        hisat2_extract_exons.py {input.genome_gtf} > {params.outdir}/genome.exon;
        hisat2_extract_splice_sites.py {input.genome_gtf} > {params.outdir}/genome.ss;
        hisat2-build -p {threads} --ss {params.outdir}/genome.ss --exon {params.outdir}/genome.exon {input.genome_fasta} {params.dbPrefix} 1>{log} 2>&1
        '''

# ----------------------------------------------------------------------------
#step2
step2_outdir = join(config["result_dir"], "align.hisat2_stringtie", \
    "step2_mapping_genome_by_hisat2")
step2_dbPrefix = join(step1_outdir, step1_dbname)
isIndexOk = join(step1_outdir, "index.ok")
if config["rna_library"] == "fr-firststrand":
    rna_library = "--rna-strandness RF"
elif config["rna_library"] == "fr-secondstrand":
    rna_library = "--rna-strandness FR"
elif config["rna_library"] == "none":
    rna_library = ""
# ----------------------------------------------------------------------------
rule align_hisat2_stringtie_step2_mapping_genome_by_hisat2:
    message:
        '''
        hisat2
        '''
    input:
        read1 = join(config["sample_dir"], "{sample}_R1.fq.gz"),
        read2 = join(config["sample_dir"], "{sample}_R2.fq.gz"),
        isIndexOk = isIndexOk,
    output:
        join(step2_outdir, "{sample}", "{sample}.sorted.bam")
    params:
        dbPrefix = step2_dbPrefix,
        rna_library = rna_library
    conda:
        py3env
    benchmark:
        join(step2_outdir, "{sample}", "benchmark.txt")
    log:
        join(step2_outdir, "{sample}", "hisat2.log")
    threads:
        5
    shell:
        '''
        hisat2 --dta {params.rna_library} -p {threads} -x {params.dbPrefix} -1 {input.read1} -2 {input.read2} 2>{log} | samtools sort -@ {threads} -o {output} 1>>{log} 2>&1
        '''

# ----------------------------------------------------------------------------
#step3
step3_inputdir = step2_outdir
step3_outdir = join(config["result_dir"], "align.hisat2_stringtie", \
                                "step3_assembly_transcript_by_stringtie")
step3_cuffPrefix = "cufcompF"
# ----------------------------------------------------------------------------
rule align_hisat2_stringtie_step3_assembly_transcript_by_stringtie:
    message:
        '''
        stringtie
        '''
    input:
        join(step3_inputdir, "{sample}", "{sample}.sorted.bam")
    output:
        join(step3_outdir, "{sample}", "{sample}.gtf")
    params:
        gff = config["genome_gtf"],
        step3_cuffPrefix = join(step3_outdir, "{sample}", step3_cuffPrefix),
        cuffout = join(step3_outdir, "{sample}", step3_cuffPrefix + ".combined.gtf")
    conda:
        py3env
    benchmark:
        join(step3_outdir, "{sample}", "benchmark.txt")
    log:
        join(step3_outdir, "{sample}", "stringtie.log")
    threads: 5
    shell:
        '''
        stringtie -p {threads} -G {params.gff} -o {output} {input} 1>{log} 2>&1;
        cuffcompare -r {params.gff} -o {params.step3_cuffPrefix} {output};
        less {params.cuffout}|perl -ne '/transcript_id "(\S+)";.*class_code "(\S)"/;print "$1\\t$2\\n"'|sort -u|cut -f2|sort|uniq -c|perl -ne 's/^\s+//;print "- $_"' >{params.step3_cuffPrefix}.classcode.stat.txt;
        '''

# ----------------------------------------------------------------------------
#step4
step4_inputdir = step3_outdir
step4_outdir = join(config["result_dir"], "align.hisat2_stringtie", "step4_merge_stringtieResult_by_stringtieMerge")
step4_cuffPrefix = "cufcompF"
# ----------------------------------------------------------------------------
rule align_hisat2_stringtie_step4_merge_stringtieResult_by_stringtieMerge:
    message:
        '''
        stringtieMerge
        '''
    input: 
        expand(join(step4_inputdir, "{sample}", "{sample}.gtf"), sample = SAMPLES),
    output:
        mergedGtf = join(step4_outdir, "merged.gtf"),
        cuffcomparedGtf = join(step4_outdir, step4_cuffPrefix + ".combined.gtf")
    params:
        gff = config["genome_gtf"],
        step4_cuffPrefix = join(step4_outdir, step4_cuffPrefix),
        cuffout = join(step3_outdir, step4_cuffPrefix + ".combined.gtf")
    conda:
        py3env
    benchmark:
        join(step4_outdir, "benchmark.txt")
    log:
        join(step4_outdir, "merged.log")
    shell:
        '''
        stringtie --merge -G {params.gff} -o {output.mergedGtf} {input} 1>{log} 2>&1;
        cuffcompare -r {params.gff} -o {params.step4_cuffPrefix} {output.mergedGtf} 1>>{log} 2>&1;
        less {output.cuffcomparedGtf}|perl -ne '/transcript_id "(\S+)";.*class_code "(\S)"/;print "$1\\t$2\\n"'|sort -u|cut -f2|sort|uniq -c|perl -ne 's/^\s+//;print "- $_"' >{params.step4_cuffPrefix}.classcode.stat.txt;
        '''
