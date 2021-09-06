
rule gff3repair:
    input:
        genome_gtf = config['genome_gtf'],
        genome_fasta = config['genome_fasta']
    output:
        genome_gff = join(config['result_dir'], 'database', 'data', 'genome.gff3'),
    conda: 'envs/agat.yaml'
    log: join(config['result_dir'], 'database', 'data', 'gff3.log'),
    shell:
        '''
        agat_convert_sp_gxf2gxf.pl -g results/result/transcript_assembly/transcript_assembly_by_stringtie/merged.gtf -o results/result/transcript_assembly/transcript_assembly_by_stringtie/merged.gff
        #python gffrename.py Pvo.gff mRNA_old2new.csv >PVDB.genome.gff
        #python fastarename.py Function_annotation.xls mRNA_old2new.csv Function_annotation.csv xingyouteng.gene.final.changeid.pep >PVDB.pep.fa
        #python fastarename.py Function_annotation.xls mRNA_old2new.csv Function_annotation.csv xingyouteng.gene.final.changeid.cds >PVDB.cds.fa
        #python gffanno.py Function_annotation.xls mRNA_old2new.csv PVDB.genome.gff >PVDB.genome.note.gff
        '''