
summary of reads produced
=========================

数据质量情况详细内容如下：

#. Sample name：样品名称
#. Raw reads：统计原始序列数据，以四行为一个单位，统计每个文件的测序序列的个数
#. Clean reads：计算方法同 Raw reads，只是统计的文件为过滤后的测序数据,后续的生物信息分析都是基于 Clean reads
#. Clean bases：测序序列的个数乘以测序序列的长度，并转化为以G为单位
#. Error rate：碱基测序错误率
#. Q20、Q30：分别表示 Phred 数值大于20、30的碱基占总体碱基的百分比
#. GC content：碱基 G 和 C 的数量总和占总的碱基数量的百分比