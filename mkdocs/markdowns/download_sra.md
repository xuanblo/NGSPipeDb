conda create -n sradownload
conda activate sradownload
mamba install sra-tools=2.10.1 -c bioconda

mamba install -y -c hcc aspera-cli=3.9.1