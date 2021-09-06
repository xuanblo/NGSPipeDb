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
