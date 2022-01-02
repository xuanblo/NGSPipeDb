# bsa analysis

!!! note "A typical flow of bsa analysis with reference is shown in the figure below"
    <figure markdown> 
        ![Dummy image](imgs/bsa.png){ width="800" }
        <figcaption>bsa pipeline</figcaption>
    </figure>

Two way to run:

=== "multi-parameter way"

    ```
    ngspipedb runpipe ngspipe-bsa-medicago -n ngspipe-bsa --resultdirname result -d ../test_pipeline --genomeFasta ../testdata_ngspipe-rnaseq-basic/genome/chr19.fa --genomeAnno ../testdata_ngspipe-rnaseq-basic/genome/GRCm38.83.chr19.gtf --samplefile ../testdata_ngspipe-rnaseq-basic/rawdata/sample.csv --rawreadsdir ../testdata_ngspipe-rnaseq-basic/rawdata --snaketype p --readsprefix _R{}.fq.gz -j 1
    ```