# Exogenous TE insertion

> - ITIS: tool to Identify TE Insertion Sites in genome
> - https://github.com/Chuan-Jiang/ITIS
> - https://pubmed.ncbi.nlm.nih.gov/25887332/

Two way to run:

=== "multi-parameter way"

    ```
    ngspipedb runpipe ngspipe-tnt-medicago -n ngspipe-tnt --resultdirname result -d ../ --genomeFasta ../genome/scf003.fa --genomeAnno ../genome/scf003.sorted.gff --samplefile ../rawdata/sample.csv --rawreadsdir ../rawdata --snaketype p --reads_prefix _R{}_paired.fastq.gz -j 20
    ```

    - `ngspipe-tnt-medicago` your project name
    - `-n ngspipe-tnt` pipeline name
    - `--genomeFasta ../genome/scf003.fa` give a genome fasta file path, see file format [fasta](https://en.wikipedia.org/wiki/FASTA_format)
    - `--genomeAnno ../genome/scf003.sorted.gff` give a genome annotaion file path [gtf](https://genome.ucsc.edu/FAQ/FAQformat.html#format4)/[gff](https://genome.ucsc.edu/FAQ/FAQformat.html#format3)
    - `--samplefile ../rawdata/sample.csv` give a sample file path, which has one column

=== "use configfile way"

    1. start a project in current directory:
        ```shell
        ngspipedb startproject medicago_tnt_zhao -n ngspipe-tnt
        ```
        - PROJECTNAME: medicago_tnt_zhao
        - pipeline: ngspipe-tnt
        - `tree medicago_tnt_zhao`
            ```
            ├── database
            ├── genome
            ├── ngsdb_config.yaml
            ├── ngspipe_config.yaml
            └── rawdata
                └── sample.csv
            3 directories, 3 files
            ```
        
    2. copy your file to directory:
        ```
        ├── database
        │   ├── exogenous
        │   │   ├── medicago_mere1.fa
        │   │   └── medicago_tnt1.fa
        ├── genome
        │   ├── GRCm38.83.chr19.gtf
        │   ├── chr19.fa
        │   └── chr19.fa.fai
        └── rawdata
            ├── control-0_R1.fq.gz
            ├── control-0_R2.fq.gz
            ├── control-1_R1.fq.gz
            ├── control-1_R2.fq.gz
            ├── control-2_R1.fq.gz
            ├── control-2_R2.fq.gz
            ├── sample.csv
            ├── treated-0_R1.fq.gz
            ├── treated-0_R2.fq.gz
            ├── treated-1_R1.fq.gz
            ├── treated-1_R2.fq.gz
            ├── treated-2_R1.fq.gz
            └── treated-2_R2.fq.gz
        ```
        - database/exogenous: the directory of exogenous sequence
        - genome: the genome file
        - rawdata: resequencing data
    
    3. modify configfile file `medicago_tnt_zhao/ngspipe_config.yaml`:

        ```yaml
        #---------------------------
        # medicago tnt1
        #---------------------------

        ## 1.reference ##
        genomeAnno_path: "genome/GRCm38.83.chr19.gtf" # gene annotation file, can be gtf or gff format
        genomeFasta_path: "genome/chr19.fa" # genome sequence, fasta format
        exogenous_seq_path: "database/exogenous/medicago_tnt1.fa" # medicago_mere1.fa or medicago_tnt1.fa

        ## 2.raw reads data ##
        sample_path: "rawdata/sample.csv" # sample file
        rawreads_dir: "rawdata" # sample file directory
        read1Suffix: "_R1.fq.gz" # fastq suffix, read1
        read2Suffix: "_R2.fq.gz"

        ## 3.output directory ##
        results_name: "results"

        ## 4.notice ##
        # if the string is 'nobody', ngspipe will not send email
        # modify 'noboby' to 'xxx@qq.com' or 'xxx@qq.com,yyy@qq.com' to send email
        email_addr: 'nobody'

        # choose where to stop your pipeline
        target: all #

        #----------------------------------
        # Configuration for sampling data
        #----------------------------------

        # for test the pipe, you can choose to the part of the input file
        # which sampling method do you want to use?
        sampling_method: links # links or head or tail or seqkit_number or seqkit_proportion
        # Default is links (ues the whole data of sample); head (use first sampling_range line in every sample),tail (use last sampling_range line in every sample); seqkit_number (number of reads); seqkit_proportion (percentage of reads)
        # and how many reads file line or reads number or reads proportion do you want to use?
        sampling_value: 80000 # for head and tail, this value is line number; for number, this value is reads number; for proportion, this value is percentage

        samples_num: all # all or interger
        # Default is all (use all samples), give a sample number, must less than real sample number, for example 6

        #----------------------------------
        # Configuration for Quality Control
        #----------------------------------
        # which qc method do you want to use?

        qc_method: trim-galore # trim-galore or trimmomatic or fastqc

        #----------------------------------
        # Configuration for Genome cut
        #----------------------------------
        # which part of genome do you want to use?
        sub_genome: all # all (whole genome) or chr1:10000-20000 (seq name:start position-end position) or chr1 (seq name)
        ```
    
    4. run command:
        ```shell
        ngspipedb runpipe medicago_tnt_zhao -n ngspipe-tnt -c medicago_tnt_zhao/ngspipe_config.yaml -j 1
        ```