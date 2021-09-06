import sys
import re
import os

replict_num = 3
proption = 0.5
#new_sufix = 'fastq.gz'

SAMPLE_NAME = re.compile(r'(\S+)_(\S+\.fq\.gz)')

for sample in sys.argv[1:]:
    matObj = SAMPLE_NAME.match(sample)
    if matObj:
        name, sufix = matObj.groups()
        #print(name, sufix)
        for i in range(replict_num):
            new_sample = name + '-' + str(i) + '_' + sufix
            command = 'gunzip -c {sample} | seqkit sample -s {seed} -p {proption} -o {new_sample}'.format(sample=sample,seed=10+i, proption=proption,new_sample=new_sample)
            print(command)
            os.system(command)