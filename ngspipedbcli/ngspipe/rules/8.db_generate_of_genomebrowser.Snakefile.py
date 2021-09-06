

rule index_genome:
    message:
        '''
        ------------------------------
        8. igv index genome
        ------------------------------
        '''
    input:
        genome = config['genomeFasta']
    output:
        genome = join(annotation_gbrowse_outdir, "genome.fa")
    log:
        join(annotation_gbrowse_outdir, "genome.log")
    shell:
        '''
        cp {input.genome} {output.genome} 2>{log};
        samtools faidx {output.genome} 2>>{log};
        '''

rule igv_annotation_gtf:
    message:
        '''
        ------------------------------
        8. igv annotation file
        ------------------------------
        '''
    input:
        refgtf = config['genomeAnno'],
    output:
        sortedGtf = join(annotation_gbrowse_outdir, 'annotation.sorted.gtf'),
        gtfgzip = join(annotation_gbrowse_outdir, 'annotation.sorted.bgzip')
    benchmark:
        join(annotation_gbrowse_outdir, "benchmark.txt")
    log:
        join(annotation_gbrowse_outdir, "run.log")
    shell:
        '''
        bedtools sort -i {input.refgtf} >{output.sortedGtf} 2>{log};
        bgzip -c {output.sortedGtf} >{output.gtfgzip} 2>>{log};
        tabix -p gff {output.gtfgzip} 2>>{log};
        '''