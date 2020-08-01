---
layout:     post   				    # 使用的布局（不需要改）
title:    A Single-Cell RNA Sequencing Profiles the Developmental Landscape of Arabidopsis Root  				# 标题 
subtitle:   My first post
date:       2020.08.01 				# 时间
author:      		Gavin				# 作者
header-img: img/post-bg-2015.jpg 	#这篇文章标题背景图片
catalog: true 						# 是否归档
tags:								#标签
    - Single-Cell
    - Plant
    - Paper
    
---
<!-- TOC -->

- [Introduction](#introduction)
  - [Structure](#structure)
  - [Mechanism](#mechanism)
- [Results](#results)
  - [scRNA-Seq and Identification of Root Cell Clusters](#scrna-seq-and-identification-of-root-cell-clusters)
    - [scRNA-Seq](#scrna-seq)
    - [Discovery](#discovery)
  - [Reconstruction of the Continuous Differentiation Trajectory of ROOT Cells](#reconstruction-of-the-continuous-differentiation-trajectory-of-root-cells)
    - [Whole root](#whole-root)
    - [Proximal Meristerm](#proximal-meristerm)
    - [Root Cap](#root-cap)
  - [Differential Hormonal Response and Ion Assimilation Pattern in Root Cell Atlas](#differential-hormonal-response-and-ion-assimilation-pattern-in-root-cell-atlas)
  - [Identification of a Role of Cytokinin Response in Lateral Root Cap Development](#identification-of-a-role-of-cytokinin-response-in-lateral-root-cap-development)
- [Discussion](#discussion)

<!-- /TOC -->
# Introduction

## Structure

![image-20200717212832441](https://i.loli.net/2020/07/17/Goe1ifTy6QDWl8Y.png)



- Root apical meristem (RAM) resides in the root tip

- The core of RAM is quiescent center (QC), known as stem cell niche (SCN), or surrounding stem cells

- other cell type

  1. epidermis
  2. hair cell
  3. root cap
  4. cortex
  5. endodermis
  6. stele
     - pericycle
     - phloem
     - xylem
     - procambium

## Mechanism

  - an intersection of high auxin concentration together with endodermal transcription factors controls the localization of the QC
  - TF WOX5, only detected in QC, maintains stem cell pluripotency by repressing the expression of CDF4, which encodes differentiation factor
  - TF PLT guides the transition from cell division to cell expansion and differentiation

    1. PLT proteins form a gradient by mitotic distribution and short distance cell-to-cell movement, with **highest level at root tip**

    2. The PLT gradient coordinates with root cell division, growth, and differentiation by **activating cell division and growth-related genes and repressing differentiation genes**
- spatial distribution of auxin and cytokinin, in coordination with PEAR and microRNA165/6-targeted HD-ZIP III transcription factors, forms a regulatory network in regulating cambium stem cell population and procambium cell development

# Results

## scRNA-Seq and Identification of Root Cell Clusters

![image-20200717220931599](https://i.loli.net/2020/07/17/uX7dJmDUTohFC9K.png)

### scRNA-Seq

- About 15 000 cells were initially mixed with 10xGenomics single-cell reaction regents for scRNA-seq assay.

- 7695 cells with 23 161 genes used for further analysis
- 24 clusters

### Discovery

1. Some cell types are scattered while other cell types were grouped together
2. A given root cell type mayconsist of different subpopulations

3. We identified several putative novel cell/subcell types in roots (clusters 6, 11, 16, 17, 22, and 24)

## Reconstruction of the Continuous Differentiation Trajectory of ROOT Cells 

### Whole root

![image-20200717222345131](https://i.loli.net/2020/07/17/BxvS3h7FZoIt4pG.png)

RC, root cap; CE, cell elongation; RH, root hair; CS, endodermal cells with casparian strip; PM, proximal meristem; M, meristem (SCN and PM).

- UMAP revealed two continuous trajectories of root
  development based on this niche: one led to the differentiation of root cap, lateral root cap, epidermis, and root hairs toward distal end and another led to the differentiation of multiple stele cell types toward the proximal end

### Proximal Meristerm

<img src="https://i.loli.net/2020/07/17/zVFpKOo2E5tyvhk.png" alt="image-20200717223419916" style="zoom: 200%;" />

SC, stem cell; PM, proximal meristem.



![image-20200717223436898](https://i.loli.net/2020/07/17/QZbtcu9fow7Y1i6.png)

### Root Cap

![image-20200717230102700](https://i.loli.net/2020/07/17/kloO8dRDKzMbgna.png)

LRC, lateral root cap; CRC, columella root cap.



- This trajectory supports the notion that LRC, epidermis, and root hair originate from the epidermal/lateral root cap initials

## Differential Hormonal Response and Ion Assimilation Pattern in Root Cell Atlas

![image-20200717231444129](https://i.loli.net/2020/07/17/9joIusg8cJaCHTZ.png)

- a single-cell root cell atlas is powerful for functional dissection of the biological processes related to hormonal response and ion assimilation at single-cell resolution.

## Identification of a Role of Cytokinin Response in Lateral Root Cap Development

![image-20200717232410231](https://i.loli.net/2020/07/17/FckMiafR5VXjOeE.png)

- arr1 arr10 arr12 triple mutant, which was largely defective in cytokinin signaling output, showed a reduced number of LRC cell layers

- these results demonstrate the functional significance of high cytokinin response during root epidermal cell differentiation and LRC development.

# Discussion

- further analyses lack crucial information on cells’ environments and locations

- the Plant Cell Atlas Project (PCA)

- The preparation of sufficient single cells for scRNA-seq is still challenging in plants.

- Due to the low number of QC cells in roots, we were unable to faithfully identify a corresponding cell cluster in our scRNA-seq dataset