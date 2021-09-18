# https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html

rule trinity_assembly1:
    '''
    Document: https://github.com/trinityrnaseq/trinityrnaseq/wiki

    '''
    input:
        read1 = expand(join(qc_outdir, "{sample}", "{sample}.cleanR1.fq.gz"), sample = SAMPLES),
        read2 = expand(join(qc_outdir, "{sample}", "{sample}.cleanR2.fq.gz"), sample = SAMPLES),
    output: 
        transcript_fa = join(transcript_assembly_outdir, 'trinity_assembly1', 'Trinity.fasta')
    params:
        trinityfolder = join(transcript_assembly_outdir, 'trinity_assembly1')
    threads: 30
    log:
        trinitysample = join(transcript_assembly_outdir, 'trinity_assembly1', 'run.log.txt')
    run:
        import os
        left_fqs = ','.join(input.read1)
        right_fqs = ','.join(input.read2)
        #print(left_fqs, right_fqs)
        command = "Trinity --seqType fq --max_memory 100G --left {left_fqs} --right {right_fqs} --output {trinityfolder} --CPU {threads} 1>{log} 2>&1;".format(left_fqs=left_fqs, right_fqs=right_fqs, trinityfolder=params.trinityfolder, threads=threads, log=log)
        print(command)
        os.system(command)

rule trinity_assembly2:
    '''
    Document: https://github.com/trinityrnaseq/trinityrnaseq/wiki

    '''
    input:
        read1 = expand(join(qc_outdir, "{sample}", "{sample}.cleanR1.fq.gz"), sample = SAMPLES),
        read2 = expand(join(qc_outdir, "{sample}", "{sample}.cleanR2.fq.gz"), sample = SAMPLES),
    output: 
        samples_file = join(transcript_assembly_outdir, 'trinity_assembly2', 'Trinity.samples.tsv'),
        transcript_fa = join(transcript_assembly_outdir, 'trinity_assembly2', 'Trinity.fasta'),
        gene_trans_map = join(transcript_assembly_outdir, 'trinity_assembly2', 'Trinity.fasta.gene_trans_map'),
    params:
        trinityfolder = join(transcript_assembly_outdir, 'trinity_assembly2'),
        qc_outdir = qc_outdir,
    threads: 30
    log:
        trinitysample = join(transcript_assembly_outdir, 'trinity_assembly2', 'run.log.txt')
    shell:
        '''
        python {snake_dir}/scripts/sample4trinity_read_condition.py {config[condition_path]} {params.qc_outdir} >{output.samples_file} 2>{log};
        Trinity --seqType fq --max_memory 100G --samples_file {output.samples_file} --output {params.trinityfolder} --CPU {threads} 1>{log} 2>&1;
        '''

# downstrean analysis

# Transcript Quantification

rule prepare_the_reference_for_alignment_and_abundance_estimation:
    input:
        transcript_fa = join(transcript_assembly_outdir, 'trinity_assembly2', 'Trinity.fasta'),
        gene_trans_map = join(transcript_assembly_outdir, 'trinity_assembly2', 'Trinity.fasta.gene_trans_map'),
    output:
        prepareok = touch(join(transcript_assembly_outdir, 'prepare_the_reference_for_alignment_and_abundance_estimation', "prepare.ok")),
    log:
        abundance_log = join(transcript_assembly_outdir, 'prepare_the_reference_for_alignment_and_abundance_estimation', 'prepare.log.txt')
    threads:
        30
    params:
        est_method = 'RSEM',
        aln_method = 'bowtie',
    shell:
        '''
        ## Just prepare the reference for alignment and abundance estimation
        align_and_estimate_abundance.pl --transcripts {input.transcript_fa} --est_method {params.est_method} --aln_method {params.aln_method} --gene_trans_map {input.gene_trans_map} --prep_reference 1>{log} 2>&1;
        '''


rule Estimating_transcript_abundance:
    input:
        transcript_fa = join(transcript_assembly_outdir, 'trinity_assembly2', 'Trinity.fasta'),
        gene_trans_map = join(transcript_assembly_outdir, 'trinity_assembly2', 'Trinity.fasta.gene_trans_map'),
        prepareok = join(transcript_assembly_outdir, 'prepare_the_reference_for_alignment_and_abundance_estimation', "prepare.ok"),
        read1 = join(qc_outdir, "{sample}", "{sample}.cleanR1.fq.gz"),
        read2 = join(qc_outdir, "{sample}", "{sample}.cleanR2.fq.gz"),
    output:
        abundance_gene = join(transcript_assembly_outdir, 'Estimating_transcript_abundance', '{sample}', "{}.isoforms.results".format('RSEM')),
    log:
        abundance_log = join(transcript_assembly_outdir, 'Estimating_transcript_abundance', '{sample}', 'abundance.log.txt')
    threads:
        10
    params:
        est_method = 'RSEM',
        aln_method = 'bowtie',
        abundance_gene_outdir = join(transcript_assembly_outdir, 'Estimating_transcript_abundance', '{sample}'),
    shell:
        '''
        ## Run the alignment and abundance estimation (assumes reference has already been prepped, errors-out if prepped reference not located.)
        align_and_estimate_abundance.pl --transcripts {input.transcript_fa} --thread_count {threads} --seqType fq --left {input.read1} --right {input.read2} --est_method {params.est_method} --aln_method {params.aln_method} --gene_trans_map {input.gene_trans_map} --output_dir {params.abundance_gene_outdir} 1>>{log} 2>&1;
        '''

rule Build_Transcript_and_Gene_Expression_Matrices:
    input:
        abundance_gene = expand(join(transcript_assembly_outdir, 'Estimating_transcript_abundance', '{sample}', "{}.isoforms.results".format('RSEM')), sample = SAMPLES),
        gene_trans_map = join(transcript_assembly_outdir, 'trinity_assembly2', 'Trinity.fasta.gene_trans_map'),
    output:
        gene_matrix_tsv = join(quantify_outdir, "trinity.gene.counts.matrix"),
        gene_matrix_csv = join(quantify_outdir, "gene.csv"),
    params:
        out_prefix = join(quantify_outdir, "trinity"),
    log:
        abundance_log = join(quantify_outdir, 'matrix.log.txt')
    shell:
        '''
        abundance_estimates_to_matrix.pl --est_method RSEM --gene_trans_map {input.gene_trans_map} --cross_sample_norm TMM --name_sample_by_basedir --out_prefix {params.out_prefix} {input.abundance_gene} 1>{log} 2>&1;
        less {output.gene_matrix_tsv}|tr '\t' ',' >{output.gene_matrix_csv} 2>>{log};
        '''

rule SuperTranscripts:
    input:
        transcript_fa = join(transcript_assembly_outdir, 'trinity_assembly2', 'Trinity.fasta'),
    output:
        supertranscript_fa = join(transcript_assembly_outdir, 'SuperTranscripts', 'trinity_genes.fasta'),
        supertranscript_gff = join(transcript_assembly_outdir, 'SuperTranscripts', 'trinity_genes.gtf'),
    params:
        out_prefix = join(transcript_assembly_outdir, 'SuperTranscripts', "trinity_genes"),
    log:
        supertranscript_log = join(transcript_assembly_outdir, 'SuperTranscripts', 'run.log.txt'),
    shell:
        '''
        Trinity_gene_splice_modeler.py --trinity_fasta {input.transcript_fa} --out_prefix {params.out_prefix} 1>{log} 2>&1;
        '''