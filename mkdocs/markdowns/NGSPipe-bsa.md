https://github.com/FenglabBioinfo/DMAP_M2

perl DMAP_M2.pl --vcf ./Example/allsnp.vcf  --BulkMut PM --BulkWild PW --win 1000000 --step 500000 --out ./Example/pp2a
解压后，里面有一个DMAP_M2.pl的脚本，就是分析用的脚本。运行命令，就是上面这句话。

allsnp.vcf 就是样本的变异信息文件。
 --BulkMut PM --BulkWild PW 对应突变体和野生型池的名字。
--win 1000000 --step 500000 计算平均值的窗口大小和滑窗步进大小

allsnp.vcf 就是样本的变异信息文件。
 --BulkMut PM --BulkWild PW 对应突变体和野生型池的名字。
--win 1000000 --step 500000 计算平均值的窗口大小和滑窗步进大小

程序只能在Linux下运行，要安装好perl、R和ggplot2包。

一般都要修改这个子函数sub filter_step2{｝
修改过滤条件，如两个池中index值大小。