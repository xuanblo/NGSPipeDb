---
title: 
date: 2021-03-07 09:51:19
tags:
  - Markdown
  - rnaseq
categories: module
identification of novel transcripts, identification of expressed, alternative splicing, and for the detection of gene fusion events.

---

# Reference-free RNA-seq analysis use NGSPipe

#### 

!!! Info inline end
    If this is your first time using NGSPipe, then we strongly recommend that you start by running test data. If you already have experience with NGSPipe, we suggest you can go straight to the custom data section.

Reference genome-free - no genome assembly for the species of interest is available. In this case one would need to assemble the reads into transcripts using de novo approaches. This type of RNAseq is as much of an art as well as science because assembly is heavily parameter-dependent and difficult to do well.
In this lesson we will focus on the Reference genome-based type of RNA seq.

---

??? note "A typical flow of transcriptome analysis with reference is shown in the figure below"
    ![img](../imgs/RNA-seq-analysis-flow-chart-An-example-RNA-seq-analysis-workflow-is-depicted-for-a_W640.jpeg)

## Quick Start - One time installation of components necessary for RNA-Seq analysis on test data <a name="QuickStarted"></a>

