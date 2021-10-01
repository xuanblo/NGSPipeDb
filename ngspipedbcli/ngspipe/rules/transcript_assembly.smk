# 4. transcript assembly
transcript_assembly_method = config['transcript_assembly_method'] # star
stringtie_rna_library = "" # --rf: Assumes a stranded library fr-firststrand. --fr: Assumes a stranded library fr-secondstrand.
transcript_assembly_outdir = join(config["resultsDir"], "transcript_assembly")

rule transcript_assembly_by_stringtie:
    message:
        '''
        stringtie was used to assemble transcript
        '''
    input:
        sorted_bam = join(junction_align_outdir, junction_align_method, 'align', "{sample}", "{sample}.sorted.bam"),
        genomeAnno = config["genomeAnno"]
    output:
        gtf = join(transcript_assembly_outdir, 'stringtie', 'assembly', "{sample}", "{sample}.gtf")
    benchmark:
        join(transcript_assembly_outdir, 'stringtie', 'assembly', "{sample}", "benchmark.txt")
    log:
        join(transcript_assembly_outdir, 'stringtie', 'assembly', "{sample}", "stringtie.log.txt")
    threads: 5
    shell:
        '''
        stringtie {stringtie_rna_library} -p {threads} -G {input.genomeAnno} -o {output.gtf} {input.sorted_bam} 1>{log} 2>&1;
        '''

# ----------------------------------------------------------------------------
# merged_gtf
# ----------------------------------------------------------------------------
rule merge_stringtieResult_by_stringtieMerge:
    message:
        '''
        stringtieMerge merge all gtf file together
        '''
    input: 
        gtfs = expand(join(transcript_assembly_outdir, 'stringtie', 'assembly', "{sample}", "{sample}.gtf"), sample = SAMPLES),
        genomeAnno = config["genomeAnno"]
    output:
        mergedGtf = join(transcript_assembly_outdir, "merged.gtf")
    benchmark:
        join(transcript_assembly_outdir, "benchmark.txt")
    log:
        join(transcript_assembly_outdir, "merged.log")
    shell:
        '''
        stringtie --merge -G {input.genomeAnno} -o {output.mergedGtf} {input.gtfs} 1>{log} 2>&1;
        '''

rule transcript_assembly:
    message:
        '''
        transcript_assembly
        '''
    input:
        mergedGtf = join(transcript_assembly_outdir, "merged.gtf")
    output:
        transcript_assembly_ok = touch(join(flag_outdir, 'transcript_assembly.ok'))