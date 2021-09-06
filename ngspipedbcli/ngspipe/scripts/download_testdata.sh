# get directory
dest_forder=$1

if [ ! -d $dest_forder ]; then
  mkdir $dest_forder
fi

echo -e "\033[32m Entering forder: \033[31m ${dest_forder} \033[0m, start to download: \033[0m"
cd $dest_forder

# Download the file of interest (here using a loop)
# Note that it can be interesting to store the STDERR of wget.
echo -e "\033[32m Downloading raw reads... \033[0m"
# --no-clobber 不要重复下载已存在的文件
for i in control_R1 control_R2 treated_R1 treated_R2;
  do [[ -e ${i}.fq.gz ]] || wget --no-clobber http://www.liu-lab.com/ngspipedb/testdata/${i}.fq.gz;
done

# Download chr19 sequence (mm10 version)
echo -e "\033[32m Downloading reference sequence... \033[0m"
[[ -e chr19.fa ]] || wget --no-clobber http://www.liu-lab.com/ngspipedb/testdata/chr19.fa.gz
[[ -e chr19.fa ]] || gunzip chr19.fa.gz # uncompress the file

echo -e "\033[32m Downloading gff file... \033[0m"
curl -C - http://www.liu-lab.com/ngspipedb/testdata/GRCm38.83.chr19.gtf.gz | gunzip -c > GRCm38.83.chr19.gtf

echo -e "\033[32m Downloading sample list file... \033[0m"
wget --no-clobber http://www.liu-lab.com/ngspipedb/testdata/samples.xls

echo -e "\033[32m Downloading condition list file... \033[0m"
wget --no-clobber http://www.liu-lab.com/ngspipedb/testdata/condition.xls

`python ../ngspipe/scripts/generate_replicat.py control_R1.fq.gz control_R2.fq.gz treated_R1.fq.gz treated_R2.fq.gz`
#`rm -f control_R1.fq.gz control_R2.fq.gz treated_R1.fq.gz treated_R2.fq.gz`

# cd back
cd -

echo -e "\033[32m Finished. Back to current direcotory \033[0m"