# 1. download ngspipedb to anywhere you want
# git clone git://github.com/xuanblo/NGSPipeDb.git && mv NGSPipeDb mouse_transcriptome_analysis && cd mouse_transcriptome_analysis

platform=$1

# 2. install conda for your platform
echo -e "\033[32m [Step1:] install conda for your platform \033[0m"

if [ $platform == 'linux' ]
then
    [[ -e /tmp/`whoami`_Miniconda3-latest-Linux-x86_64.sh ]] || wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/`whoami`_Miniconda3-latest-Linux-x86_64.sh
    bash /tmp/`whoami`_Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda3
elif [ $platform == 'macos' ]
then
    [[ -e /tmp/`whoami`_Miniconda3-latest-Linux-x86_64.sh ]] || wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O /tmp/`whoami`_Miniconda3-latest-Linux-x86_64.sh bash /tmp/`whoami`_Miniconda3-latest-Linux-x86_64.sh -b -f -p ~/miniconda3
else
    exit -1
fi

# conda init
echo -e "\033[32m [Step2:] conda init and update \033[0m"
~/miniconda3/bin/conda init bash
source ~/miniconda3/etc/profile.d/conda.sh && conda update conda -y && conda install mamba -c conda-forge -y

# 3. create conda environment # 如果已经存在的话会询问
echo -e "\033[32m [Step3:] Create conda env \033[0m"
conda create -c conda-forge -c bioconda --name ngspipe-rnaseq snakemake=5.30.2 python=3.8 seqkit=0.14.0 -y 

# 4. update ensential packages
echo -e "\033[32m [Step4:] update ensential packages \033[0m"
conda env update -n ngspipe-rnaseq --file ngspipe/envs/requirements_rnaseq.yaml --prune

# 5. enter conda env
echo -e "\033[32m [Step5:] activate conda env \033[0m"
conda activate ngspipe-rnaseq

# 6. download testfile
echo -e "\033[32m [Step6:] download test data \033[0m"
bash ngspipe/scripts/download_testdata.sh testdata

# 7. run snakemake
echo -e "\033[32m [Step7:] run RNA-Seq analysis \033[0m"
snakemake -s ngspipe/1.1.rnaseq_analysis_reference_basic.Snakefile.py --configfile ngspipe/config/rnaseq.config.yaml --dag|dot -Tpng > dag.png
snakemake -s ngspipe/1.1.rnaseq_analysis_reference_basic.Snakefile.py --configfile ngspipe/config/rnaseq.config.yaml -p -j 10

# 8. generate report
#snakemake --snakefile ngspipe/1.1.rnaseq_analysis_reference_basic.Snakefile.py --configfile ngspipe/config/rnaseq.config.yaml --report results/report/report.html

# 9. exit env
#conda deactivate

# rm -rf ~/miniconda
# tac ~/.bashrc|sed "1,15{d}"|tac > /tmp/bashrc && mv /tmp/bashrc ~/.bashrc
# sed -i "$(($(wc -l < ~/.bashrc) - 15)),\$d" ~/.bashrc