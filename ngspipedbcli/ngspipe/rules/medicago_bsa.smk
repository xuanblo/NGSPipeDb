bsa_outdir = join(config["resultsDir"], "bsa")

# QTL-BSA: A Bulked Segregant Analysis and Visualization Pipeline for QTL-seq

rule download_itis_code:
    output:
        code_dir = directory(join(bsa_outdir, 'software'))
    shell:
        '''
        mkdir -p {output.code_dir};
        cd {output.code_dir};
        git clone git://github.com/FenglabBioinfo/DMAP_M2.git;
        '''

rule genome_size:
    input:
        genomeFasta_out = join(sub_genome_dir, 'sub_genome.fa'),
    output:
        genome_size = join(bsa_outdir, 'index', 'genome_size.txt')
    log: join(bsa_outdir, 'index', 'run.log.txt')
    shell:
        '''
        faSize {input.genomeFasta_out} >{output.genome_size} 2>{log};
        '''

rule vcf_merge:
    input:
        vcf = expand(join(variation_outdir, "{sample}", "{sample}.vcf"), sample=SAMPLES)
    output:
        vcf_all_samples = join(bsa_outdir, "vcf_merge", "all.vcf"),
        vcf_all_samples_rename = join(bsa_outdir, "vcf_merge", "all_rename.vcf"),
        genome_size = config['genomeFasta'] + '.fai'
    log:
    run:
        import os
        vcfs = input.vcf
        vcf_params = ['-V ' + i for i in vcfs[0:2]]
        vcf_params_str = ' '.join(vcf_params)
        command1 = "samtools faidx {}".format(config['genomeFasta'])
        command2 = "picard CreateSequenceDictionary R={}".format(config['genomeFasta'])
        command3 = 'gatk3 -T CombineVariants {} -o {} -R {}'.format(vcf_params_str, output.vcf_all_samples, config['genomeFasta'])
        os.system(command1)
        os.system(command2)
        os.system(command3)

        new_vcf = open(output.vcf_all_samples_rename, 'w')
        with open(output.vcf_all_samples, 'r') as f:
            for line in f.readlines():
                if line.startswith('#CHROM'):
                    items = line.split('\t')
                    new_items = items[0:9] + SAMPLES[0:2]
                    new_line = '\t'.join(new_items) + '\n'
                    new_vcf.write(new_line)
                else:
                    new_vcf.write(line)


rule bsa_find_by_DMAP:
    message:
        '''
        https://github.com/FenglabBioinfo/DMAP_M2
        '''
    input:
        vcf = join(bsa_outdir, "vcf_merge", "all_rename.vcf"),
        genome_size = config['genomeFasta'] + '.fai',
        code_dir = join(bsa_outdir, 'software'),
    output:
        variation = join(bsa_outdir, 'DMAP', 'variation_info.txt'),
    threads: 30
    params:
        output_dir = join(bsa_outdir, 'DMAP')
    log: join(bsa_outdir, 'DMAP', 'run.log.txt')
    shell:
        '''
        perl {input.code_dir}/DMAP_M2/DMAP_M2.pl --vcf {input.vcf} --index {input.genome_size} --BulkMut {SAMPLES[0]} --BulkWild {SAMPLES[1]} --win 1000000 --step 500000 --out {params.output_dir} 1>{log} 2>&1;
        '''

rule bsa_merge:
    input:
        variation = join(bsa_outdir, 'DMAP', 'variation_info.txt'),
    output:
        bsa_merge_ok = touch(join(flag_outdir, 'bsa_merge.ok'))