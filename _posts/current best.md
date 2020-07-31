---
layout:     post   				    # 使用的布局（不需要改）
title:     Current best practices in single-cell RNA-seq analysis a tutorial  				# 标题 
subtitle:   #My first post
date:       2020.8.1 				# 时间
author:      						# 作者
header-img: img/post-bg-2015.jpg 	#这篇文章标题背景图片
catalog: true 						# 是否归档
tags:								#标签
    - Single-Cell
    - Bioinformatics
    - Paper
    
---

<!-- TOC -->

- [Introduction](#introduction)
    - [The challenges to standardization](#the-challenges-to-standardization)
    - [Key elements of an experimental scRNA-seq workflow](#key-elements-of-an-experimental-scrna-seq-workflow)
        - [single-cell isolation](#single-cell-isolation)
            - [methods](#methods)
            - [error](#error)
        - [library construction](#library-construction)
        - [sequencing](#sequencing)
            - [Without UMI](#without-umi)
            - [With UMI](#with-umi)
- [Content](#content)
    - [Pre-processing and visualization](#pre-processing-and-visualization)
        - [Tools](#tools)
        - [Procedure](#procedure)
            - [read quality control (QC)](#read-quality-control-qc)
            - [demultiplexing](#demultiplexing)
            - [genome alignment](#genome-alignment)
            - [quantification](#quantification)
    - [Quality control](#quality-control)
        - [Goal](#goal)
        - [Procedure](#procedure-1)
            - [Checking  the integrity of cells](#checking--the-integrity-of-cells)
                - [Tools](#tools-1)
            - [Checking the level of transcripts](#checking-the-level-of-transcripts)
            - [Correcting the ambient gene expression](#correcting-the-ambient-gene-expression)
                - [Tools](#tools-2)
        - [Pitfalls & recommendations](#pitfalls--recommendations)
    - [Normalization](#normalization)
        - [Tools](#tools-3)
            - [Cellular count data](#cellular-count-data)
                - [Count depth scaling (CPM)](#count-depth-scaling-cpm)
                - [Extension of CPM](#extension-of-cpm)
                - [Non-linear normalization methods](#non-linear-normalization-methods)
                    - [Tools](#tools-4)
                - [Appropriate normalization method for a specific dataset](#appropriate-normalization-method-for-a-specific-dataset)
                    - [Full-length](#full-length)
                    - [3' enrichment](#3-enrichment)
            - [Gene counts](#gene-counts)
            - [log(x+1) transformation](#logx1-transformation)
                - [advantage](#advantage)
                - [disadvantage](#disadvantage)
        - [Pitfalls & recommendations:](#pitfalls--recommendations)
    - [Data correction and integration](#data-correction-and-integration)
        - [Goal](#goal-1)
        - [Methods](#methods)
            - [Regressing out biological effects](#regressing-out-biological-effects)
                - [Remove the effects of the cell cycle on the transcriptome](#remove-the-effects-of-the-cell-cycle-on-the-transcriptome)
                - [Aspects considered](#aspects-considered)
            - [Regressing out technical effects](#regressing-out-technical-effects)
            - [Batch effects and data intergration](#batch-effects-and-data-intergration)

<!-- /TOC -->
# Introduction

## The challenges to standardization

- Growing number of
  analysis methods and exploding
  dataset sizes
- Analysis tools for scRNA-seq data are written in a variety of
  programming languages

Instead of targeting a standardized analysis pipeline, we outline current best practices
and common tools independent of programming language 



##  Key elements of an experimental scRNA-seq workflow

*Ziegenhain et al (2017);Macosko et al (2015); Svensson et al (2017).*

![image-20200710145425486](https://i.loli.net/2020/07/10/eZjoE6BamJcfQq1.png)

<center>figure 1</center>



### single-cell isolation

#### methods

1. plate-based techniques isolate cells into wells on a plate
2. droplet-based methods rely on capturing each cell in its own microfluidic droplet

#### error

errors can occur that lead to
multiple cells being captured together (doublets or multiplets), nonviable
cells being captured, or no cell being captured at all (empty
droplets/wells).

### library construction

(Library construction is the process in which the intracellular mRNA is captured, reverse-transcribed to cDNA molecules and  amplified)

1. Each well or droplet contains the necessary chemicals to break down the cell membranes and perform library construction

2. the mRNA from each cell can
   be labelled with a well- or droplet-specific cellular barcode

3. experimental protocols also label captured molecules (mRNA) with a unique molecular identifier (UMI).

   (**UMIs allow us to distinguish between amplified copies of the same mRNAmolecule and reads from separate mRNA molecules transcribed from the same gene.**)


### sequencing

#### Without UMI

Sequencing produces **read data**, which undergo quality control, grouping based on their assigned barcodes (demultiplexing) and alignment in read processing pipelines.

#### With UMI

For UMI-based protocols, read data can be further demultiplexed to produce counts of captured mRNA molecules
(**count data**).

# Content

## Pre-processing and visualization

### Tools

Cell Ranger (Zheng et al, 2017), indrops (Klein et al, 2015),
SEQC (Azizi et al, 2018), or zUMIs (Parekh et al, 2018)

### Procedure

#### read quality control (QC)

#### demultiplexing

（assigning reads to their cellular barcodes and mRNA molecules of origin ）

#### genome alignment

#### quantification

## Quality control

![image-20200710150515627](https://i.loli.net/2020/07/10/GolXmDpF4YyqibN.png)

<center>figure2</center>

### Goal

 Ensuring that all cellular barcode data correspond to **viable** cells.

### Procedure

#### Checking  the integrity of cells

There QC covariates（jointly）

1. the number of counts per barcode (count depth)
2. the number of genes per barcode
3. the fraction of counts from mitochondrial genes per barcode

##### Tools

*DoubletDecon: preprint: DePasquale et al, 2018; Scrublet: Wolock et al, 2019;Doublet Finder: McGinnis et al, 2018*

#### Checking the level of transcripts

<font color=red>A guideline to setting this threshold is to use the minimum cell cluster size that is of interest and leaving some leeway for dropout effects</font>

#### Correcting the ambient gene expression

**Ambient gene expression** refers to counts that do not originate from a barcoded cell, but from other lysed cells whose mRNA contaminated the cell suspension prior to library construction.

##### Tools

*SoupX (preprint: Young & Behjati, 2018)*

*Angelidis et al, 2019*

### Pitfalls & recommendations

-  Perform QC by finding outlier peaks in the number of genes, the count depth and the fraction of mitochondrial reads. Consider these covariates jointly instead of separately.
-  Be as permissive of QC thresholding as possible, and revisit QC if downstream clustering cannot be interpreted.
-  If the distribution of QC covariates differ between samples, QC thresholds should be determined separately for each sample to account for sample quality differences as in Plasschaert et al (2018).

## Normalization

### Tools

- bulk gene expression

  preprint: Pachter, 2011; Dillies et al, 2013

- single-cell data

  Lun et al, 2016a; Vallejos et al, 2017

  ### Procedure

####  Cellular count data

##### Count depth scaling (CPM)

  - also referred to as “counts per million” or CPM normalization.
  - This protocol comes from bulk expression analysis and normalizes count data using a so-called **size factor** proportional to the count depth per cell.
  - assuming that all cells in the dataset initially contained an equal number of mRNA molecules and count depth differences arise only due to sampling.

##### Extension of CPM

- high-count filtering CPM

  *Weinreb et al (2018)*

- *Scran’s pooling-based size factor estimation method (Lun et al, 2016a)*

  size factors are estimated based on a **linear regression** over genes, after cells are pooled to avoid technical dropout effects.

##### Non-linear normalization methods

- particularly relevant for **plate-based** scRNA-seq data, which tend to have batch effects between plates.

- Such an approach can combine technical and biological data correction (e.g. batch correction or correction for cell cycle effects) with count depth normalization.

###### Tools

- Mayer et al (2018) 

  it a negative binomial model to count data, using technical covariates such as the read depth and the number of counts per gene to fit the model parameters.The residuals of the model fit serve as a normalized quantification of gene expression.

- *Cole et al, 2019*

  outperform global scaling methods especially in situations with strong batch effects 

##### Appropriate normalization method for a specific dataset

- *Cole et al (2019)*

  their **scone** tool should be used to select the appropriate normalization method for a specific dataset

- ScRNA-seq techniques can be divided into **full-length** and **3' enrichment** methods (Svensson et al, 2017; Ziegenhain et al, 2017).

###### Full-length

- benefit from normalization methods that take into account gene length (e.g. Patel et al, 2014; Kowalczyk et al, 2015; Soneson & Robinson, 2018)

- TPM normalization (Li et al, 2009)

###### 3' enrichment

scran for normalization of non-full-length datasets

#### Gene counts

- Gene normalization constitutes scaling gene counts to have zero mean and unit variance (z scores)

- <font color=blue>There is currently no
  consensus on whether or not to perform normalization over genes.</font>

1. pro:

Seurat tutorials (Butler et al, 2018) generally
apply gene scaling

2. con 

the authors of the Slingshot method opt against
scaling over genes in their tutorial (Street et al, 2018).

- The preference between the two choices revolves around whether all genes should be weighted equally for downstream analysis, or whether the magnitude of expression of a gene is an informative proxy for the importance of the gene.
- opt to refrain from scaling over genes in this tutorial.

#### log(x+1) transformation

##### advantage

- effect

1. distances between log-transformed expression values represent log fold changes, which are the canonical way to measure changes in expression.
2. log transformation mitigates (but does not remove) the mean–variance relationship in single-cell data (Brennecke et al, 2013)
3. log transformation reduces the skewness of the data to approximate the assumption of many downstream analysis tools that the data are **normally distributed.**

- This usefulness is highlighted by down stream applications for **differential expression testing** (Finak et al, 2015; Ritchie et al, 2015) or **batch correction** (Johnson et al, 2006; Buttner et al, 2019) that use

##### disadvantage

introduce spurious differential expression effects into the
data (preprint: Lun, 2018). 

This effect is particularly pronouncedwhen normalization size factor distributions differ strongly between tested groups.

### Pitfalls & recommendations:



- We recommend scran for normalization of non-full-length datasets. An alternative is to evaluate normalization approaches via scone especially for plate-based datasets. 
- Full-length scRNA-seq protocols can be corrected for gene length using bulk methods.
-  Normalized data should be log(x+1)-transformed for use with downstream analysis methods that assume data are normally distributed.

## Data correction and integration

### Goal

- Data correction targets further technica and biological covariates such as **batch, dropout, or cell cycle effects.**

- the decision of which covariates to consider will depend on the intended downstream analysis.

### Methods

#### Regressing out biological effects

##### Remove the effects of the cell cycle on the transcriptome

- tools

implemented in the Scanpy and Seurat platforms (Butler
et al, 2018; Wolf et al, 2018) 

 in specialized packages with more complex mixture models such as scLVM (Buettner et al, 2015) or fscLVM (Buettner et al, 2017).

- also be used to regress out other known biological
  effects such as mitochondrial gene expression, which is interpreted as an indication of cell stress.

##### Aspects considered 

- cell cycle signals can also be informative of the biology

- correcting for one process may unintentionally mask the signal of another.
- correcting for cell size via normalization, or dedicated tools such as cgCorrect (Blasi et al, 2017), also partially corrects for cell cycle effects in scRNA-seq data.

#### Regressing out technical effects

- Most prominent technical covariates are **count depth and batch**
- Regressing out count depth effects can improve the performance of trajectory inference algorithms, which rely on finding transitions between cells

- An alternative to regression-based strategies for removing count effects is to use a more rigorous normalization procedure such as downsampling or non-linear normalization methods

#### Batch effects and data intergration

- Three levels of batch effects

1. between groups of cells in an experiment,

**batch correction**

2.between experiments performed in the same laboratory or between datasets from different laboratories.

**data integration.**

